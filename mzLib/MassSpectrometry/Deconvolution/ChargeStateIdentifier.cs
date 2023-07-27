#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;
using Chemistry;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MathNet.Numerics.Statistics;


namespace MassSpectrometry;

public class ChargeStateIdentifier : ClassicDeconvolutionAlgorithm
{
    private ChargeStateDeconvolutionParams DeconvolutionParams { get; }
    public ChargeStateIdentifier(DeconvolutionParameters deconParameters) : base(deconParameters)
    {
        DeconvolutionParams = (ChargeStateDeconvolutionParams)deconParameters;
    }
    public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrumToDeconvolute, MzRange range)
    {
        return DeconvolutePrivateFast(spectrumToDeconvolute,
            range, DeconvolutionParams.EnvelopeThreshold)
            .OrderByDescending(i => i.Score)
            .DistinctBy(i => new
            {
                i.MonoisotopicMass,
                i.Charge, 
                i.Score
            })
            .ToList();
    }

    
    /// <summary>
    /// The function that actually does the work in deconvolution, optimized for speed. 
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="deconvolutionRange"></param>
    /// <param name="spectralSimMatchThresh"></param>
    /// <returns></returns>
    internal IEnumerable<IsotopicEnvelope> DeconvolutePrivateFast(MzSpectrum scan, MzRange deconvolutionRange, double spectralSimMatchThresh)
    {
        // I think the isotopic envelope has to be a List<IsotopicEnvelope>. 
        // And charge states will get added to the isotopic envelope object. 
        // 
        ConcurrentDictionary<double, List<IsotopicEnvelope>> ieHashSet = new(new DoubleEqualityComparer());
        double[] masses = GenerateMassesIndex(DeconvolutionParams.MinimumMassDa, DeconvolutionParams.MaximumMassDa,
            DeconvolutionParams.PeakMatchPpmTolerance); 

        var output = PreFilterMzVals(scan.XArray,
            scan.YArray, DeconvolutionParams.MinCharge, 
            DeconvolutionParams.MaxCharge,
            DeconvolutionParams.MinimumMassDa,
            DeconvolutionParams.MaximumMassDa,
            masses
            );

        var ladder = CreateChargeStateLadders(output, DeconvolutionParams.MinCharge, DeconvolutionParams.MaxCharge,
            scan.FirstX!.Value, scan.LastX!.Value);

        Parallel.ForEach(ladder, (m) =>
        {
            var index = MatchChargeStateLadder(scan, m);
            var ladderMatch = TransformToChargeStateLadderMatch(index, scan, m, DeconvolutionParams.PeakMatchPpmTolerance);

            var successfulMatch = ScoreChargeStateLadderMatch(ladderMatch, scan);

            if (successfulMatch)
            {
                FindIsotopicEnvelopes(ladderMatch!, scan, deconvolutionRange, ieHashSet, DeconvolutionParams.EnvelopeThreshold);
            }
        });

        return ieHashSet.Values.SelectMany(i => i);
    }
    
    /// <summary>
    /// Given a charge state ladder match and the original data, adds 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    /// <param name="range"></param>
    /// <param name="ieHashSet"></param>
    /// <param name="minimumThreshold"></param>
    internal void FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan, MzRange range, 
        ConcurrentDictionary<double, List<IsotopicEnvelope>> ieHashSet, double minimumThreshold)
    {
        //TODO: Needs to be refactored to remove side effects. 
        this.spectrum = scan;
        List<double> chargesList = match.PeakList.Select(i => i.Charge).ToList(); 
        for (int i = 0; i < chargesList.Count; i++)
        {
            int charge = (int)Math.Round(chargesList[i]);
            if (range.Contains(chargesList[i]))
            {
                double[] neutralXArray = scan.XArray.Select(j => j.ToMass(charge)).ToArray();

                MzSpectrum neutralMassSpectrum = new(neutralXArray, scan.YArray, true);

                double maxMassToTake = allMasses[match.MassIndex].Max();
                double minMassToTake = allMasses[match.MassIndex].Min();

                MzRange newRange = new(minMassToTake, maxMassToTake);

                var envelope = FillIsotopicEnvelopeByBounds(match, neutralMassSpectrum, newRange, charge);
                if (envelope != null)
                {
                    // need to perform a scoring if the envelope only consists of low resolution data. 
                    RescoreIsotopicEnvelope(envelope);
                    if (envelope.Score >= minimumThreshold)
                    {
                        if (ieHashSet.ContainsKey(match.MonoisotopicMass))
                        {
                            if (ieHashSet.TryGetValue(match.MonoisotopicMass, out var tempList))
                            {
                                tempList.Add(envelope);
                            }
                        }
                        else
                        {
                            ieHashSet.TryAdd(match.MonoisotopicMass, new List<IsotopicEnvelope> { envelope });
                        }
                    }
                }
            }
        }
    }
    /// <summary>
    /// Generates an sequence of values where the n+1 value is a the 2 * ppm tolerance away from the previous value. 
    /// </summary>
    /// <param name="minMass"></param>
    /// <param name="maxMass"></param>
    /// <param name="ppmTolerance"></param>
    /// <returns></returns>
    internal static double[] GenerateMassesIndex(double minMass, double maxMass, double ppmTolerance = 15)
    {
        List<double> masses = new(); 
        // put mass in the center of the ppm error range. 

        double value = minMass;
        double epsilon = ppmTolerance / 1e6; 
        while (value < maxMass)
        {
            double newValue = value * (epsilon + 1)/(1 - epsilon);
            masses.Add(newValue);
            value = newValue;
        }
        return masses.ToArray();
    }

    /// <summary>
    /// Calculates theoretical charge state ladders. 
    /// </summary>
    /// <param name="indexOfMaxIntensityPeak">Index derived from the original, full scan. </param>
    /// <param name="mzValueArray"></param>
    /// <param name="minCharge"></param>
    /// <param name="maxCharge"></param>
    /// <param name="minMzValAllowed"></param>
    /// <param name="maxMzValAllowed"></param>
    /// <returns></returns>
    internal static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(int indexOfMaxIntensityPeak,
        double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed)
    {
        double mzVal = mzValueArray[indexOfMaxIntensityPeak];
        for (int i = minCharge; i <= maxCharge; i++)
        {
            double tempMass = mzVal.ToMass(i);
            List<(int charge, double mz)> tempLadder = new();
            for (int j = maxCharge; j > minCharge; j--)
            {
                double mz = tempMass.ToMz(j);
                if (mz <= maxMzValAllowed && mz >= minMzValAllowed)
                {
                    tempLadder.Add((j, tempMass.ToMz(j)));
                }
            }
            // clean up the ladder ;
            yield return new ChargeStateLadder(tempMass, tempLadder.Select(k => k.mz).ToArray());
        }
    }

    internal static IEnumerable<double> PreFilterMzVals(double[] mzVals, double[] intensityArray, int minCharge, 
        int maxCharge, double minMass, double maxMass, double[] masses, double ppmMatchTolerance = 15)
    {
        ConcurrentDictionary<double, int> concurrentHashSet = new(new DoubleEqualityComparer());
        ConcurrentDictionary<int, double> cumulativeIntensityDict = new();
        for (int i = 0; i < masses.Length; i++)
        {
            cumulativeIntensityDict.TryAdd(i, 0d); 
        }
        

        // contains the neutral mass and the index of the neutral mass
        
        
        Parallel.For(0, mzVals.Length, j =>
        {
            for (int i = maxCharge; i >= minCharge; i--)
            {
                var testMass = mzVals[j].ToMass(i);
                if (testMass > maxMass || testMass < minMass)
                {
                    continue;
                }

                int index = GetBucket(masses, testMass);
                concurrentHashSet.TryAdd(testMass, index);
                while (true)
                {
                    if (cumulativeIntensityDict.TryUpdate(index, cumulativeIntensityDict[index] + intensityArray[j],
                            cumulativeIntensityDict[index])) break;
                }
            }
        });

        // Define the threshold value as mean + 1.5 * standard deviation of the intensity values
        var cumulativeIntensityArray = cumulativeIntensityDict.Values.ToArray();
        var meanVariance = ArrayStatistics.QuantileInplace(cumulativeIntensityArray.Where(z => z > 0).ToArray(), 0.5); 
        double threshold = meanVariance;

        //double medianCounts = ArrayStatistics.Mean(cumulativeIntensities.Where(i => i >= 1.1).ToArray());
        // return indexes where counts are bigger than the median. 
        HashSet<int> indexToKeepList = new();
        
        for (int i = 0; i < cumulativeIntensityArray.Length; i++)
        {
            if (cumulativeIntensityArray[i] >= threshold)
            {
                indexToKeepList.Add(i); 
            }
        }

        return concurrentHashSet
            .Where(i => indexToKeepList.Contains((i.Value)))
            .GroupBy(x => x.Value, new DoubleEqualityComparer())
            .ToDictionary(t => t.Key, 
                t => t.Select(r => r.Key).Average())
            .Select(i => i.Value);
    }


    internal static int GetBucket(double[] array, double value)
    {
        int index = Array.BinarySearch(array, value);
        if (index < 0)
        {
            index = ~index;
        }

        if (index >= array.Length)
        {
            index = array.Length - 1;
        }

        if (index != 0 && array[index] - value >
            value - array[index - 1])
        {
            index--;
        }
        return index;
    }
    internal static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(IEnumerable<double> neutralMass, 
        int minCharge, int maxCharge, double minMzAllowed, double maxMzAllowed)
    {
        foreach (var mass in neutralMass)
        {
            List<(int charge, double mz)> tempLadder = new();
            for (int j = maxCharge; j > minCharge; j--)
            {
                double mz = mass.ToMz(j);
                if (mz <= maxMzAllowed && mz >= minMzAllowed)
                {
                    tempLadder.Add((j, mass.ToMz(j)));
                }
            }
            yield return new ChargeStateLadder(mass, tempLadder.Select(k => k.mz).ToArray());
        }
    }

    internal static List<int> MatchChargeStateLadder(MzSpectrum scan, ChargeStateLadder ladder)
    {
        ConcurrentBag<int> output = new();
        foreach (var t in ladder.MzVals)
        {
            int index = GetBucket(scan.XArray, t);
            output.Add(index);
        }

        return output.ToList();
    }

    internal static ChargeStateLadderMatch TransformToChargeStateLadderMatch(List<int> ladderToIndexMap, 
        MzSpectrum scan, ChargeStateLadder ladder, double ppmMatchTolerance)
    {
        List<double> listMzVals = new();
        List<double> listIntVals = new();
        for (int i = 0; i < ladderToIndexMap.Count; i++)
        {
            listMzVals.Add(scan.XArray[ladderToIndexMap[i]]);
            listIntVals.Add(scan.YArray[ladderToIndexMap[i]]);
        }
        List<double> chargesOfMatchingPeaks = listMzVals.Select(i => ladder.Mass / i).ToList();

        return new ChargeStateLadderMatch(listMzVals, listIntVals, 
            chargesOfMatchingPeaks, ladder, ppmMatchTolerance); 
    }
    


    /// <summary>
    /// Performs spectral similarity calculation between the theoretical isotopic envelopes and the experimental data.
    /// </summary>
    /// <param name="envelope"></param>
    internal void RescoreIsotopicEnvelope(IsotopicEnvelope envelope)
    {
        if (envelope == null) return;

        double diff = envelope.MonoisotopicMass + diffToMonoisotopic[envelope.MassIndex] - allMasses[envelope.MassIndex][0]; 
        double[] theoreticalMzs = allMasses[envelope.MassIndex].Select(i => (i + diff).ToMz(envelope.Charge)).ToArray();
        double maxIntensity = allIntensities[envelope.MassIndex].Max(); 
        double[] theoreticalIntensities = allIntensities[envelope.MassIndex]
            .Select(i => i / maxIntensity)
            .ToArray();

        var theoreticalSpectrum = new MzSpectrum(theoreticalMzs, theoreticalIntensities, true)
            .FilterByY(0.001, 1.0);
        var spectrum0 = new MzSpectrum(theoreticalSpectrum.Select(i => (i.Mz)).ToArray(),
            theoreticalSpectrum.Select(i => i.Intensity).ToArray(), true); 

        MzSpectrum experimentalSpectrum = new MzSpectrum(envelope.Peaks.Select(i => i.mz.ToMz(envelope.Charge)).ToArray(),
            envelope.Peaks.Select(i => i.intensity).ToArray(), true);

        SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, spectrum0,
            SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, DeconvolutionParams.PeakMatchPpmTolerance,
            true);

        double? score = similarity.SpectralContrastAngle();
        if (score.HasValue)
        {
            envelope.Rescore(score.Value);
        }
        else
        {
            envelope.Rescore(0);
        }
    }
    
    internal bool ScoreChargeStateLadderMatch(ChargeStateLadderMatch match, MzSpectrum scan)
    {
        GetMassIndex(match, scan);
        
        match.CompareTheoreticalNumberChargeStatesVsActual();
        if (match.PercentageMzValsMatched < DeconvolutionParams.PercentageMatchedThresh) return false;

        match.CalculateChargeStateScore();
        if (match.SequentialChargeStateScore < DeconvolutionParams.SequentialChargeStateDiff) return false;
        
        match.CalculateEnvelopeScore();
        if (match.EnvelopeScore < DeconvolutionParams.EnvelopeScoreThresh) return false;

        CalculateMonoMass(match);

        return true; 
    }
    
    internal void CalculateMonoMass(ChargeStateLadderMatch match)
    {
        match.MonoGuesses = match.PeakList
            .Select(i => i.NeutralMass - diffToMonoisotopic[match.MassIndex])
            .ToList();
    }

    /// <summary>
    /// Adds the mass index to the ChargeStateLadderMatch object. Must be run from the deconvolution to have access to the precalculated averagine isotopic envelopes. 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    public void GetMassIndex(ChargeStateLadderMatch match, MzSpectrum scan)
    {
        // binary search to get the mass index and assign it to match. 
        this.spectrum = scan;
        // get the index of the averagine model. 
        int massIndex = Array.BinarySearch(mostIntenseMasses, match.TheoreticalLadder.Mass);
        if (massIndex < 0)
        {
            massIndex = ~massIndex;
        }

        if (massIndex >= mostIntenseMasses.Length)
        {
            massIndex = mostIntenseMasses.Length - 1;
        }

        if (massIndex != 0 && mostIntenseMasses[massIndex] - match.TheoreticalLadder.Mass >
            match.TheoreticalLadder.Mass - mostIntenseMasses[massIndex - 1])
        {
            massIndex--;
        }
        match.MassIndex = massIndex;
    }
    /// <summary>
    /// get the range of the theoretical isotopic envelopes and then pulls all the peaks within that range from the original data. 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    /// <param name="isolationRange"></param>
    /// <param name="chargeState"></param>
    /// <returns></returns>
    public IsotopicEnvelope FillIsotopicEnvelopeByBounds(ChargeStateLadderMatch match, MzSpectrum scan, MzRange isolationRange, int chargeState)
    {
        List<(double, double)> listOfPeaks = scan.Extract(isolationRange).Select(i => (i.Mz, i.Intensity)).ToList();
        double totalIntensity = listOfPeaks.Sum(i => i.Item2);
        if (listOfPeaks.Any())
        {
            return new(listOfPeaks, match.MonoisotopicMass, chargeState,
                totalIntensity, match.StdDev, match.MassIndex);
        }

        return null; 
    }
}

public class DoubleEqualityComparer : IEqualityComparer<double>
{
    public bool Equals(double a, double b)
    {
        return Math.Round(a, 2) == Math.Round(b, 2); 
    }

    public int GetHashCode(double value)
    {
        return Math.Round(value, 2).GetHashCode();
    }
}