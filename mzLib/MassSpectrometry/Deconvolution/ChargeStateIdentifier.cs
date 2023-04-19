#nullable enable
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;
using System.Runtime.ExceptionServices;
using System.Runtime.InteropServices;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;

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


// link to where I got this code in stackOverFlow: https://stackoverflow.com/a/18950824
    public static IEnumerable<Tuple<int, double>> FindLocalMaxima(IEnumerable<double> array, int windowSize)
    {
        // Round up to nearest odd value
        windowSize = windowSize - windowSize % 2 + 1;
        int halfWindow = windowSize / 2;

        int index = 0;
        var before = new Queue<double>(Enumerable.Repeat(double.NegativeInfinity, halfWindow));
        var after = new Queue<double>(array.Take(halfWindow + 1));

        foreach (double d in array.Skip(halfWindow + 1).Concat(Enumerable.Repeat(double.NegativeInfinity, halfWindow + 1)))
        {
            double curVal = after.Dequeue();
            if (before.All(x => curVal > x) && after.All(x => curVal >= x))
            {
                yield return Tuple.Create(index, curVal);
            }

            before.Dequeue();
            before.Enqueue(curVal);
            after.Enqueue(d);
            index++;
        }
    }

    public IEnumerable<Tuple<int, double>> ConvertToIndexedTuple(IEnumerable<double> array)
    {
        int index = 0; 
        foreach (double d in array)
        {
            yield return Tuple.Create(index, d);
            index++; 
        }
    }

    private ChargeStateDeconvolutionParams DeconvolutionParams { get; set; }
    public ChargeStateIdentifier(ChargeStateDeconvolutionParams deconParameters) : base(deconParameters)
    {
        DeconvolutionParams = deconParameters;
    }
    // int indexOfMaxIntensityPeak,
    // double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed,
    // double adductMass = Constants.ProtonMass
    public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
    {
        HashSet<double> seenMzValues = new();

        // run until the length of the usedMzValues does not change 
        int usedMzValuesLength1 = 0; 
        int usedMzValuesLength2 = 1;
        int continuationSum = 1;
        int iterations = 0;

        // peaks in order of intensity

        List<(int index, double mz, double intensity)> mzIntensityPairs =
            (from i in Enumerable.Range(0, spectrum.XArray.Length)
                select (i, spectrum.XArray[i], spectrum.YArray[i]))
            .OrderByDescending(i => i.Item3)
            .ToList();

        // go to next if it is found in seen 
        int indexer = 0;
        while (indexer < spectrum.XArray.Length && iterations < 150)
        {
            iterations++;
            // check if peak has been found in seen and is above minimum threshold of intensity
            if (mzIntensityPairs[indexer].intensity / mzIntensityPairs.Max(i => i.intensity) < 0.005 
                || seenMzValues.Contains(mzIntensityPairs[indexer].mz))
            {
                indexer++; 
                if(indexer == spectrum.XArray.Length) break;
                
                continue; 
            }

            var chargeStateLadders = CreateChargeStateLadders(
                mzIntensityPairs[indexer].Item1, 
                spectrum.XArray,
                DeconvolutionParams.MinCharge,
                DeconvolutionParams.MaxCharge, spectrum.XArray[0],
                spectrum.XArray[^1]).ToList();
            
            MatchChargeStateLadders(spectrum, chargeStateLadders.ToList(), 
                DeconvolutionParams.PeakMatchPpmTolerance, out List<List<int>> ladderToIndicesMaps);

            var chargeStateLadderMatches = TransformToChargeStateLadderMatch(ladderToIndicesMaps, 
                mzIntensityPairs.OrderBy(i => i.Item1)
                    .Select(i => (i.Item2, i.Item3)).ToList(), chargeStateLadders).ToList();
            ChargeStateLadderMatch? bestScoringChargeStateLadderMatch = ScoreChargeStateLadderMatches(chargeStateLadderMatches, spectrum);
            if (bestScoringChargeStateLadderMatch == null)
            {
                continue; 
            }

            var isotopicEnvelopesIntermediate = FindIsotopicEnvelopes(bestScoringChargeStateLadderMatch!, spectrum).ToList();

            foreach (var envelope in isotopicEnvelopesIntermediate)
            {
                if (envelope == null) continue; 
                if (envelope.Peaks.Count < DeconvolutionParams.MinCharge - 1) continue;
                foreach (var mz in envelope.Peaks.Select(i => i.mz))
                {
                    seenMzValues.Add(mz);
                }
                yield return envelope;
            }

            indexer++; 
        }
    }
    public ChargeStateLadderMatch? ScoreChargeStateLadderMatches(List<ChargeStateLadderMatch> ladderMatches, MzSpectrum scan)
    {

        var orderByIntensity = ladderMatches
            .Where(i => i.MatchingMzPeaks.Count > i.TheoreticalLadder.MzVals.Length / 3)
            .OrderByDescending(i => i.IntensitiesOfMatchingPeaks.Sum())
            .ToList();
        // check for harmonics
        double[] ladderMasses = orderByIntensity.Select(i => i.TheoreticalLadder.Mass).ToArray();

        List<(int, int)> indicesList = new List<(int, int)>(); 

        for (int i = 0; i < ladderMasses.Length; i++)
        {
            for (int j = 0; j < ladderMasses.Length; j++)
            {
                // i will always be the larger value. 
                double harmonic = ladderMasses[i] / ladderMasses[j];
                if (Math.Abs(harmonic - 2) <= 0.05 || Math.Abs(harmonic-3) <= 0.05 || Math.Abs(harmonic - 4) <= 0.05)
                {
                    indicesList.Add((i, j));
                }
            }
        }

        if (indicesList.Any())
        {
            indicesList.Select(i => i.Item2)
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

    public double PenalizeNonConsecutiveChargeStates(double[] chargeStates)
    {
        if(chargeStates.Length == 0) return 0;
        double[] diffs = new double[chargeStates.Length - 1];
        for (int i = 0; i < diffs.Length; i++)
        {
            diffs[i] = chargeStates[i] - chargeStates[i + 1];
        }

        double expectedSum = chargeStates.Length;
        double actualSum = diffs.Sum();
        return expectedSum / actualSum;
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
        return match.IntensitiesOfMatchingPeaks.Where(i => i >= threshold).Sum();
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
                match.TheoreticalLadder.Mass, (int)Math.Round(match.ChargesOfMatchingPeaks[i]), 5.0,
                3.0, monoisotopicMassPrediction); 

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