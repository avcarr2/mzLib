﻿#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;
using Chemistry;
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
        return DeconvolutePrivateFast(spectrumToDeconvolute, range, DeconvolutionParams.EnvelopeThreshold);
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
        ConcurrentDictionary<double, IsotopicEnvelope> ieHashSet = new(new DoubleEqualityComparer());

        // slow step about 200 ms
        var output = PreFilterMzVals(scan.XArray,
            scan.YArray, DeconvolutionParams.MinCharge, DeconvolutionParams.MaxCharge,
            DeconvolutionParams.MinimumMassDa,
            DeconvolutionParams.MaximumMassDa,
            DeconvolutionParams.PeakMatchPpmTolerance,
            DeconvolutionParams.DeltaMass);

        var ladder = CreateChargeStateLadders(output, DeconvolutionParams.MinCharge, DeconvolutionParams.MaxCharge,
            scan.FirstX!.Value, scan.LastX!.Value);

        Parallel.ForEach(ladder, (m) =>
        {
            var index = MatchChargeStateLadder(scan, m,
                DeconvolutionParams.PeakMatchPpmTolerance);
            var ladderMatch = TransformToChargeStateLadderMatch(index, scan, m);

            var successfulMatch = ScoreChargeStateLadderMatch(ladderMatch, scan);

            if (successfulMatch)
            {
                FindIsotopicEnvelopes(ladderMatch!, scan, deconvolutionRange, ieHashSet, DeconvolutionParams.EnvelopeThreshold);
            }
        });

        return ieHashSet.Values;
    }

    /// <summary>
    /// Provides access to the diff to monoisotopic mass value stored in the deconvolution algorithm abstract class.
    /// Used for plotting. Probably able to get rid of this at some point. 
    /// </summary>
    /// <param name="massIndex"></param>
    /// <returns></returns>
    public static double GetDiffToMonoisotopic(int massIndex)
    {
        return diffToMonoisotopic[massIndex];
    }

    public static double GetClosestBinMonoMass(int massIndex)
    {
        return allMasses[massIndex][0] - diffToMonoisotopic[massIndex];
    }
    
    internal void FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan, MzRange range, 
        ConcurrentDictionary<double, IsotopicEnvelope> ieHashSet, double minimumThreshold)
    {
        this.spectrum = scan;
        for (int i = 0; i < match.ChargesOfMatchingPeaks.Count; i++)
        {
            if (!ieHashSet.ContainsKey(match.MonoisotopicMass))
            {
                int charge = (int)Math.Round(match.ChargesOfMatchingPeaks[i]);
                if (range.Contains(match.MatchingMzPeaks[i]))
                {
                    double[] neutralXArray = scan.XArray.Select(j => j.ToMass(charge)).ToArray();

                    MzSpectrum neutralMassSpectrum = new(neutralXArray, scan.YArray, true);

                    double maxMassToTake = allMasses[match.MassIndex].Max();
                    double minMassToTake = allMasses[match.MassIndex].Min();

                    MzRange newRange = new(minMassToTake, maxMassToTake);

                    var envelope = FillIsotopicEnvelopeByBounds(match, neutralMassSpectrum, newRange, charge);
                    if (envelope != null)
                    {
                        RescoreIsotopicEnvelope(envelope);
                        if (envelope.Score >= minimumThreshold)
                        {
                            ieHashSet.TryAdd(match.MonoisotopicMass, envelope);
                        }
                    }
                }
            }
        }
        
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

    internal static IEnumerable<double> PreFilterMzVals(double[] mzVals, double[] intensityArray, int minCharge, int maxCharge, 
        double minMass, double maxMass, double ppmMatchTolerance = 15, double delta = 0.50)
    {
        double[] masses = new double[(int)((maxMass - minMass) / delta)];
        ConcurrentDictionary<double, int> concurrentHashSet = new(new DoubleEqualityComparer());
        ConcurrentDictionary<int, double> cumulativeIntensityDict = new();

        for (int i = 0; i < masses.Length; i++)
        {
            masses[i] = minMass + delta * i;
            cumulativeIntensityDict.TryAdd(i, 0);
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
        var meanVariance = ArrayStatistics.Mean(cumulativeIntensityArray.Where(z => z > 0).ToArray()); 
        double threshold = 1.5 * meanVariance;

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
            .Where(i => indexToKeepList.Contains(i.Value))
            .Select(i => i.Key);
    }
    private static int GetBucket(double[] array, double value)
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

    internal static List<int> MatchChargeStateLadder(MzSpectrum scan, ChargeStateLadder ladder, 
        double ppmMatchTolerance)
    {
        ConcurrentBag<int> output = new();
        foreach (var t in ladder.MzVals)
        {
            double tolerance = ppmMatchTolerance / 1e6 * t;
            double upperTol = t + tolerance;
            double lowerTol = t - tolerance;
            // get indices of the range of double values that contain the values in scan.Xarray 
            int upperIndex = Array.BinarySearch(scan.XArray, upperTol);
            int lowerIndex = Array.BinarySearch(scan.XArray, lowerTol);

            upperIndex = upperIndex < 0 ? ~upperIndex : upperIndex;
            lowerIndex = lowerIndex < 0 ? ~lowerIndex : lowerIndex;

            upperIndex = upperIndex >= scan.XArray.Length ? scan.XArray.Length - 1 : upperIndex;
            lowerIndex = lowerIndex == 0 ? 0 : lowerIndex;
            // search over a subset of the Xarray to speed up the attempted peak matching. 
            for (int k = lowerIndex; k <= upperIndex; k++)
            {
                if (scan.XArray[k] >= lowerTol && scan.XArray[k] <= upperTol)
                {
                    output.Add(k);
                }
            }
        }

        return output.ToList();
    }

    internal static ChargeStateLadderMatch TransformToChargeStateLadderMatch(List<int> ladderToIndexMap, 
        MzSpectrum scan, ChargeStateLadder ladder)
    {
        List<double> listMzVals = new();
        List<double> listIntVals = new();
        for (int i = 0; i < ladderToIndexMap.Count; i++)
        {
            listMzVals.Add(scan.XArray[ladderToIndexMap[i]]);
            listIntVals.Add(scan.YArray[ladderToIndexMap[i]]);
        }
        List<double> chargesOfMatchingPeaks = listMzVals.Select(i => ladder.Mass / i).ToList();

        return new ChargeStateLadderMatch()
        {
            MatchingMzPeaks = listMzVals,
            IntensitiesOfMatchingPeaks = listIntVals,
            TheoreticalLadder = ladder,
            ChargesOfMatchingPeaks = chargesOfMatchingPeaks
        };
    }
    


    /// <summary>
    /// Performs spectral similarity calculation between the theoretical isotopic envelopes and the experimental data.
    /// </summary>
    /// <param name="envelope"></param>
    internal void RescoreIsotopicEnvelope(IsotopicEnvelope envelope)
    {
        if (envelope == null) return;

        double diff = envelope.MonoisotopicMass + diffToMonoisotopic[envelope.MassIndex] - allMasses[envelope.MassIndex][0]; 
        double[] theoreticalMzs = allMasses[envelope.MassIndex].Select(i => i + diff).ToArray();
        double maxIntensity = allIntensities[envelope.MassIndex].Max(); 
        double[] theoreticalIntensities = allIntensities[envelope.MassIndex]
            .Select(i => i / maxIntensity)
            .ToArray();

        var theoreticalSpectrum = new MzSpectrum(theoreticalMzs, theoreticalIntensities, true)
            .FilterByY(0.001, 1.0);
        var spectrum0 = new MzSpectrum(theoreticalSpectrum.Select(i => i.Mz).ToArray(),
            theoreticalSpectrum.Select(i => i.Intensity).ToArray(), true); 

        MzSpectrum experimentalSpectrum = new MzSpectrum(envelope.Peaks.Select(i => i.mz).ToArray(),
            envelope.Peaks.Select(i => i.intensity).ToArray(), true);

        SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, spectrum0,
            SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, 5,
            false);

        double? score = similarity.CosineSimilarity();
        if (score.HasValue)
        {
            envelope.Rescore(score.Value);
        }
        else
        {
            envelope.Rescore(0);
        }
    }

    internal static double GetTheoreticalMostAbundantMass(int massIndex)
    {
        return allMasses[massIndex][0];
    }
    internal static double GetTheoreticalMostIntenseMassWithDiff(int massIndex, double diff)
    {
        return allMasses[massIndex][0] + diff;
    }


    internal bool ScoreChargeStateLadderMatch(ChargeStateLadderMatch match, MzSpectrum scan)
    {
        GetMassIndex(match, scan);
        List<double> monoGuesses = new();
        List<double> errorInMostIntense = new();
        //assume that the selected m / z is the most intense peak in the envelope. 

        if (match.TheoreticalLadder.Mass <= DeconvolutionParams.MinimumMassDa) return false;

        match.CompareTheoreticalNumberChargeStatesVsActual();
        if (match.PercentageMzValsMatched < 0.2) return false;

        match.CalculateChargeStateScore();
        if (match.SequentialChargeStateScore < -1.1) return false;
        match.CalculateEnvelopeScore();
        if (match.EnvelopeScore < 0.1) return false;

        for (int i = 0; i < match.MatchingMzPeaks.Count; i++)
        {
            double deChargedMass = match.MatchingMzPeaks[i]
                .ToMass((int)Math.Round(match.ChargesOfMatchingPeaks[i]));
            // get the theoretical highest peak in isotopic envelope. 
            double monoGuess = deChargedMass - diffToMonoisotopic[match.MassIndex];

            monoGuesses.Add(monoGuess);
        }

        match.MonoErrors = errorInMostIntense;
        match.MonoGuesses = monoGuesses;
        // use median because the center of the edges of the charge state distributions are likely to be 
        // more error prone. 
        match.MonoisotopicMass = monoGuesses.Median();
     
        match.ScoreByIntensityExplained(scan);
        return true;
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
                totalIntensity, 0d, match.MassIndex);
        }

        return null; 
    }
}

public class DoubleEqualityComparer : IEqualityComparer<double>
{
    public bool Equals(double a, double b)
    {
        return Math.Round(a, 3) == Math.Round(b, 3); 
    }

    public int GetHashCode(double value)
    {
        return Math.Round(value, 3).GetHashCode();
    }
}

//public static class IsotopicEnvelopeExtensions
//{
//    public static IDictionary<double, IEnumerable<IsotopicEnvelope>> ObserveAdjacentChargeStates(this IEnumerable<IsotopicEnvelope> listEnvelopes)
//    {
//        HashSet<double> massHashSet = new();
//        var orderedMonoMasses = listEnvelopes.OrderBy(i => i.MonoisotopicMass)
//            .Select(i => i.MonoisotopicMass)
//            .ToArray();
//        var orderedArrayOfEnvelopes = listEnvelopes
//            .OrderBy(i => i.MonoisotopicMass)
//            .ToList();

//        ConcurrentDictionary<double, IEnumerable<IsotopicEnvelope>> dictChargeStateEnvelopes = new();

//        int indexer = 0;
//        List<IsotopicEnvelope> tempIsotopicEnvelope = new();
//        while (indexer < orderedMonoMasses.Length)
//        {
            
//            if (massHashSet.Contains(orderedMonoMasses[indexer]))
//            {
//                tempIsotopicEnvelope.Add(orderedArrayOfEnvelopes.ElementAt(indexer)); 
//                indexer++; 
//            }
//            else
//            {
//                if (indexer == 0)
//                {
//                    massHashSet.Add(orderedMonoMasses[indexer]);
//                    indexer++; 
//                    continue;
//                }

//                dictChargeStateEnvelopes.TryAdd(orderedMonoMasses[indexer - 1], tempIsotopicEnvelope.DistinctBy(i => i.Charge));
//                tempIsotopicEnvelope = new(); 

//                massHashSet.Add(orderedMonoMasses[indexer]);
//                tempIsotopicEnvelope.Add(orderedArrayOfEnvelopes.ElementAt(indexer));
//                indexer++; 
//            }

//        }
//        return dictChargeStateEnvelopes;
//    }
//}






