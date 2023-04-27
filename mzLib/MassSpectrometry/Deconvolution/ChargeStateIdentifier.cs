#nullable enable
using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using System.Threading.Tasks;
using Accord.Statistics.Distributions.Univariate;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;
using Accord.Statistics.Testing;
using MathNet.Numerics.Statistics;

namespace MassSpectrometry;

public class ChargeStateIdentifier : ClassicDeconvolutionAlgorithm
{
    public IEnumerable<ChargeStateLadder> CreateChargeStateLadders(List<Tuple<int,double>> indexIntensityTupleList, double[] mzValueArray, int minCharge, int maxCharge, 
        double minMzValAllowed, double maxMzValAllowed, double adductMass = 1.007276)
    {
        // find the maximum intensity peak and use the index to get the mz values from the mzValueArray
        // (should be the first because its ordered by descending)
        int indexOfMostIntenseY = indexIntensityTupleList.First().Item1;
        double mzVal = mzValueArray[indexOfMostIntenseY];

        for (int i = minCharge; i <= maxCharge; i++)
        {
            double tempMass = (mzVal - adductMass) * (double)i; 
            List<(int charge, double mz)> tempLadder = new();
            for (int j = maxCharge; j > 0; j--)
            {
                tempLadder.Add((j,tempMass / (double)j + adductMass));
            }
            // clean up the ladder 
            tempLadder = tempLadder.Where(k => k.mz <= maxMzValAllowed && k.mz >= minMzValAllowed).ToList(); 
            yield return new ChargeStateLadder(tempMass, tempLadder.Select(k => k.charge).ToArray(), 
                i, tempLadder.Select(k => k.mz).ToArray());
        }
    }

    public IEnumerable<ChargeStateLadder> CreateChargeStateLadders(int indexOfMaxIntensityPeak,
        double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed, 
        double adductMass = Constants.ProtonMass)
    {
        double mzVal = mzValueArray[indexOfMaxIntensityPeak];
        for (int i = minCharge; i <= maxCharge; i++)
        {
            double tempMass = (mzVal - adductMass) * (double)i;
            List<(int charge, double mz)> tempLadder = new();
            for (int j = maxCharge; j > 0; j--)
            {
                tempLadder.Add((j, tempMass / (double)j + adductMass));
            }
            // clean up the ladder 
            tempLadder = tempLadder.Where(k => k.mz <= maxMzValAllowed && k.mz >= minMzValAllowed).ToList();
            yield return new ChargeStateLadder(tempMass, tempLadder.Select(k => k.charge).ToArray(),
                i, tempLadder.Select(k => k.mz).ToArray());
        }
    }

    public static List<(double mz,double intensity)> CreateMzValueArray(MzSpectrum spectrum, List<Tuple<int, double>> indexIntensityTupleList)
    {
        var orderedByIndex = indexIntensityTupleList.OrderBy(i => i.Item1);
        List<(double, double)> outputList = new List<(double, double)>();  
        foreach (var index in orderedByIndex)
        {
            outputList.Add((spectrum.XArray[index.Item1], index.Item2));
        }
        return outputList;
    }

    public List<(double,double)> MatchChargeStateLadders(MzSpectrum spectrum, List<ChargeStateLadder> chargeStateLadders,
        List<Tuple<int, double>> indexIntensityTupleList, double ppmMatchTolerance, out List<List<int>> ladderToIndicesMap)
    {
        var mzIntensityPairs = CreateMzValueArray(spectrum, indexIntensityTupleList).OrderBy(i => i.mz).ToList();
        ladderToIndicesMap = new List<List<int>>(); 
        foreach (var ladder in chargeStateLadders)
        {
            List<int> outputArray = new();
            for (int i = 0; i < ladder.MzVals.Length; i++)
            {
                for (int j = 0; j < mzIntensityPairs.Count; j++)
                {
                    double tolerance = ppmMatchTolerance / 1e6 * ladder.MzVals[i];
                    double upperTol = ladder.MzVals[i] + tolerance;
                    double lowerTol = ladder.MzVals[i] - tolerance;
                    if (mzIntensityPairs[j].mz >= lowerTol && mzIntensityPairs[j].mz <= upperTol)
                    {
                        outputArray.Add(j);
                    }
                }
            }
            
            ladderToIndicesMap.Add(outputArray);
        }

        return mzIntensityPairs; 
    }

    public void MatchChargeStateLadders(MzSpectrum spectrum, List<ChargeStateLadder> chargeStateLadders, 
        double ppmMatchTolerance, out List<List<int>> ladderToIndicesMaps)
    {
        ladderToIndicesMaps = new List<List<int>>();
        foreach (var ladder in chargeStateLadders)
        {
            List<int> outputArray = new();
            for (int i = 0; i < ladder.MzVals.Length; i++)
            {
                for (int j = 0; j < spectrum.XArray.Length; j++)
                {
                    double tolerance = ppmMatchTolerance / 1e6 * ladder.MzVals[i];
                    double upperTol = ladder.MzVals[i] + tolerance;
                    double lowerTol = ladder.MzVals[i] - tolerance;
                    if (spectrum.XArray[j] >= lowerTol && spectrum.XArray[j] <= upperTol)
                    {
                        outputArray.Add(j);
                    }
                }
            }
            ladderToIndicesMaps.Add(outputArray);
        }
    }
    public void MatchChargeStateLaddersFast(MzSpectrum scan, List<ChargeStateLadder> chargeStateLadders,
        double ppmMatchTolerance, out List<List<int>> ladderToIndicesMaps)
    {
        

        var ladderDict = new ConcurrentDictionary<int, List<int>>();

        Parallel.ForEach(Enumerable.Range(0, chargeStateLadders.Count), ladder =>
        {

            List<int> outputArray = new();
            for (int j = 0; j < chargeStateLadders[ladder].MzVals.Length; j++)
            {
                double tolerance = ppmMatchTolerance / 1e6 * chargeStateLadders[ladder].MzVals[j];
                double upperTol = chargeStateLadders[ladder].MzVals[j] + tolerance;
                double lowerTol = chargeStateLadders[ladder].MzVals[j] - tolerance;
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
                        outputArray.Add(k);
                    }
                }
            }

            ladderDict.TryAdd(ladder, outputArray);
        });
        ladderToIndicesMaps = ladderDict.OrderBy(i => i.Key)
            .Select(i => i.Value)
            .ToList(); 
    }

    public IEnumerable<ChargeStateLadderMatch> TransformToChargeStateLadderMatch(List<List<int>> ladderToIndicesMap, 
        List<(double,double)> mzIntensityList, List<ChargeStateLadder> ladders)
    {
        for (var index = 0; index < ladderToIndicesMap.Count; index++)
        {
            var indicesSet = ladderToIndicesMap[index];
            var ladder = ladders[index];

            List<double> listMzVals = new();
            List<double> listIntVals = new();
            foreach (var indices in indicesSet)
            {
                listMzVals.Add(mzIntensityList[indices].Item1);
                listIntVals.Add(mzIntensityList[indices].Item2);
            }

            List<double> chargesOfMatchingPeaks = listMzVals.Select(i => (ladder.Mass - Constants.ProtonMass) / i).ToList();

            ChargeStateLadderMatch ladderMatch = new()
            {
                MatchingMzPeaks = listMzVals,
                IntensitiesOfMatchingPeaks = listIntVals,
                TheoreticalLadder = ladder, 
                ChargesOfMatchingPeaks = chargesOfMatchingPeaks
            }; 
            yield return ladderMatch;
        }
    }
    

    public static void RemoveIdentifiedMzVals()
    {
        throw new NotImplementedException();
    }

    private ChargeStateDeconvolutionParams DeconvolutionParams { get; set; }
    public ChargeStateIdentifier(DeconvolutionParameters deconParameters) : base(deconParameters)
    {
        DeconvolutionParams = (ChargeStateDeconvolutionParams)deconParameters;
    }
    // int indexOfMaxIntensityPeak,
    // double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed,
    // double adductMass = Constants.ProtonMass
    
    public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum scan, MzRange range)
    {
        HashSet<double> seenMzValues = new();
        HashSet<double> forbiddenMasses = new(); 
        int iterations = 0;

        List<(int index, double mz, double intensity)> mzIntensityPairs =
            (from i in Enumerable.Range(0, scan.XArray.Length)
                select (i, scan.XArray[i], scan.YArray[i]))
            .OrderByDescending(i => i.Item3)
            .ToList();
        List<(int Index, double mz, double intensity)> filteredPairs =
            (from i in Enumerable.Range(0, scan.XArray.Length)
             where scan.XArray[i] >= range.Minimum && scan.XArray[i] <= range.Maximum
             select (i, scan.XArray[i], scan.YArray[i]))
            .OrderByDescending(i => i.Item3)
            .ToList(); 

     // go to next if it is found in seen 
     int indexer = 0;
        while (indexer < filteredPairs.Count && iterations < 500)
        {
            
            // check if peak has been found in seen and is above minimum threshold of intensity
            if (filteredPairs[indexer].intensity / filteredPairs.Max(i => i.intensity) < 0.05
                || seenMzValues.Contains(filteredPairs[indexer].mz)
                )
            {
                indexer++; 
                if(indexer == scan.XArray.Length) break;

                continue; 
            }
            
            var chargeStateLadders = CreateChargeStateLadders(
                filteredPairs[indexer].Item1, 
                scan.XArray,
                DeconvolutionParams.MinCharge,
                DeconvolutionParams.MaxCharge, scan.XArray[0],
                scan.XArray[^1]).ToList();
            
            MatchChargeStateLaddersFast(scan, chargeStateLadders.ToList(), 
                DeconvolutionParams.PeakMatchPpmTolerance, out List<List<int>> ladderToIndicesMaps);

            var chargeStateLadderMatches = TransformToChargeStateLadderMatch(ladderToIndicesMaps, 
                mzIntensityPairs.OrderBy(i => i.index)
                    .Select(i => (i.mz, i.intensity)).ToList(), chargeStateLadders).ToList();
            ChargeStateLadderMatch? bestScoringChargeStateLadderMatch = ScoreChargeStateLadderMatches(chargeStateLadderMatches, scan, 
                forbiddenMasses);

            if (bestScoringChargeStateLadderMatch == null)
            {
                iterations++;
                continue; 
            }
            forbiddenMasses.Add(bestScoringChargeStateLadderMatch.TheoreticalLadder.Mass);

            var isotopicEnvelopesIntermediate = FindIsotopicEnvelopes(bestScoringChargeStateLadderMatch!, scan).ToList();
            double? maxScore = double.MinValue; 
            if (isotopicEnvelopesIntermediate.All(i => i != null))
            {
                maxScore = isotopicEnvelopesIntermediate.Where(i => i != null)
                    .MaxBy(i => i.Score).Score;
            }

            foreach (var envelope in isotopicEnvelopesIntermediate)
            {
                if (envelope == null) continue;
                if (envelope.Peaks.Count < DeconvolutionParams.MinCharge - 1) continue;
                if (envelope.Score == maxScore)
                {
                    //seenMzValues.AddRange(envelope.Peaks.Select(i => i.mz)); 
                }
                //foreach (var mz in envelope.Peaks.Select(i => i.mz))
                //{
                //    seenMzValues.Add(mz);
                //}
                if (range == null)
                {
                    yield return envelope; 
                }
                else
                {
                    if (envelope.Peaks.Any(i => i.mz <= range.Maximum)
                        && envelope.Peaks.Any(i => i.mz >= range.Minimum))
                    {
                        yield return envelope;
                    }
                }
            }
            indexer++;
            iterations++;
        }
    }

    public void ScoreShapiroWilk()
    {

    }
    
    
    public ChargeStateLadderMatch? ScoreChargeStateLadderMatches(List<ChargeStateLadderMatch> ladderMatches, MzSpectrum scan, 
        HashSet<double> forbiddenMasses)
    {

        var orderByIntensity = ladderMatches
            //.Where(i => i.MatchingMzPeaks.Count > i.TheoreticalLadder.MzVals.Length / 3)
            .OrderByDescending(i => i.IntensitiesOfMatchingPeaks.Sum())
            .ToList();
        // check for harmonics
        double[] ladderMasses = orderByIntensity.Select(i => i.TheoreticalLadder.Mass).ToArray();

        List<(int, int)> indicesList = new List<(int, int)>();


        for (int i = 0; i < ladderMasses.Length; i++)
        {
            // generate harmonics 
            for (int j = 0; j < ladderMasses.Length; j++)
            {
                // We need to eliminate masses that have already been checked.
                // Therefore, we check to see if the ladder mass is contained by the hashset of 
                // forbidden masses and remove the index of the forbidden masses if they are present. 
                if (forbiddenMasses.Contains(ladderMasses[j]))
                {
                    indicesList.Add((i, j));
                    continue;
                }
                // i will always be the larger value. 
                double harmonic = ladderMasses[i] / ladderMasses[j];
                for (int k = this.DeconvolutionParams.MinCharge; k <= this.DeconvolutionParams.MaxCharge; k++)
                {
                    if(Math.Abs(harmonic - k) <= 0.05)
                    {
                        indicesList.Add((i, j));
                    }
                }
            }
        }

        if (indicesList.Any())
        {
            indicesList.Select(i => i.Item1)
                .Distinct()
                .OrderByDescending(i => i)
                .ForEach(i =>
                {
                    orderByIntensity.RemoveAt(i);
                });
            return orderByIntensity.OrderByDescending(i =>
            {
                double score = ScoreByIntensityExplained(i, scan) * CompareTheoreticalNumberChargeStatesVsActual(i);
                return score;
            }).First(); 
        }
        else
        {
            return orderByIntensity.OrderByDescending(i =>
            {
                var results = CalcRSquared(i);
                double slopeResults = Math.Abs(1d + results.slope);

                return slopeResults;
            }).FirstOrDefault(); 
        }
    }

    public void EstimateNormalDistributionParameters(double[] yarray, out double mean, out double variance)
    {
        double oneOverN = 1d / yarray.Length;

        double mu_intermediate = 0; 
        for (int i = 0; i < yarray.Length; i++)
        {
            mu_intermediate += yarray[i]; 
        }
        mean = oneOverN * mu_intermediate;

        double meanCopy = mean; 
        variance = yarray.Sum(i => Math.Pow(i - meanCopy, 2)) / yarray.Length;
    }

    public static double GetTheoreticalMostAbundantIsotopicPeak(int massIndex)
    {
        return mostIntenseMasses[massIndex];
    }

    public static double GetDiffToMonoisotopic(int massIndex)
    {
        return diffToMonoisotopic[massIndex];
    }

    public (double rsquared, double slope) CalcRSquared(ChargeStateLadderMatch cslm)
    {
        if (cslm.ChargesOfMatchingPeaks.Count < this.DeconvolutionParams.MinCharge) return (0, 0); 
        (double a, double b) lineParams = Fit.Line(cslm.ChargesOfMatchingPeaks.ToArray(), 
            Enumerable.Range(0, cslm.ChargesOfMatchingPeaks.Count)
                .Select(i => (double)i)
                .ToArray());

        var modeledValues = Enumerable.Range(0, cslm.ChargesOfMatchingPeaks.Count)
            .Select(i => (double)i * lineParams.b + lineParams.a).ToArray(); 
        
        return (GoodnessOfFit.CoefficientOfDetermination(modeledValues, cslm.ChargesOfMatchingPeaks), lineParams.b); 
    }

    public double CompareTheoreticalNumberChargeStatesVsActual(ChargeStateLadderMatch ladderMatch)
    {
        return (double)ladderMatch.ChargesOfMatchingPeaks.Count / (double)ladderMatch.TheoreticalLadder.MzVals.Length;
    }
    public double ScoreByIntensityExplained(ChargeStateLadderMatch match, MzSpectrum spectrum, double threshold = 0.01)
    {
        double medianSpectraValue = spectrum.YArray.Median(); 
        return match.IntensitiesOfMatchingPeaks.Where(i => i - medianSpectraValue >= threshold * spectrum.YArray.Max()).Sum();
    }


    public IEnumerable<IsotopicEnvelope> FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan)
    {
        this.spectrum = scan; 
        // get the index of the averagine model. 
        int massIndex = Array.BinarySearch(mostIntenseMasses, match.TheoreticalLadder.Mass);
        if (massIndex < 0)
        {
            massIndex = ~massIndex;
        }

        if (massIndex == mostIntenseMasses.Length)
        {
            yield return null;
        }

        if (massIndex >= mostIntenseMasses.Length)
        {
            yield break;
        }

        if (massIndex != 0 && mostIntenseMasses[massIndex] - match.TheoreticalLadder.Mass >
            match.TheoreticalLadder.Mass - mostIntenseMasses[massIndex - 1])
        {
            massIndex--;
        }

        List<double> monoisotopicMassPrediction = new(); 
        
        for (int i = 0; i < match.ChargesOfMatchingPeaks.Count; i++)
        {
            yield return FindIsotopicEnvelope(massIndex, match.MatchingMzPeaks[i], match.IntensitiesOfMatchingPeaks[i],
                match.TheoreticalLadder.Mass, (int)Math.Round(match.ChargesOfMatchingPeaks[i]), 
                5.0, 3.0, monoisotopicMassPrediction); 

        }

    }
}

public record struct ChargeStateLadder
{
    public double Mass;
    public int[] ChargeStatesContained;
    public int ChargeStateOfMostIntenseMass;
    public double[] MzVals;

    public ChargeStateLadder(double mass, int[] chargeStatesContained,
        int chargeStateOfMostIntenseMass, double[] mzVals)
    {
        Mass = mass; 
        ChargeStatesContained = chargeStatesContained;
        ChargeStateOfMostIntenseMass = chargeStateOfMostIntenseMass;
        MzVals = mzVals;
    }
}

public class ChargeStateLadderMatch
{
    public ChargeStateLadder TheoreticalLadder { get; set; }
    public List<double> MatchingMzPeaks { get; set; }
    public List<double> IntensitiesOfMatchingPeaks { get; set; }
    public List<double> ChargesOfMatchingPeaks { get; set; }

    public (double mzMostIntense, double maxIntensityPeak) GetMostIntenseTuple()
    {

        int indexOfMostIntense = IntensitiesOfMatchingPeaks.IndexOf(IntensitiesOfMatchingPeaks.Max()); 
        return (MatchingMzPeaks[indexOfMostIntense], IntensitiesOfMatchingPeaks[indexOfMostIntense]);
    }
}

public class ChargeStateDeconvolutionParams : DeconvolutionParameters
{
    public int MinCharge { get; }
    public int MaxCharge { get; }
    public double PeakMatchPpmTolerance { get; }

    public ChargeStateDeconvolutionParams(int minCharge, int maxCharge, double peakMatchTolerancePpm) : base()
    {
        MinCharge = minCharge; 
        MaxCharge = maxCharge;
        PeakMatchPpmTolerance = peakMatchTolerancePpm;
    }
}