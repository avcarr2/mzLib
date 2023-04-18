using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Easy.Common.Interfaces;

namespace MassSpectrometry;

public class ChargeStateIdentifier
{
    public static MzSpectrum DownsampleData(MzSpectrum spectrum, int downsampleInt)
    {
        return new MzSpectrum(
            spectrum.XArray.Where((_, index) => index % downsampleInt == 0)
                .Select(i => i)
                .ToArray(),
            spectrum.YArray.Where((_, index) => index % downsampleInt == 0)
                .Select(i => i)
                .ToArray(),
            true
        ); 
    }

    public static double[] MovingAverage(MzSpectrum spectrum, int window)
    {
        double[] buffer = new double[window];
        double[] output = new double[spectrum.YArray.Length]; 
        int currentIndex = 0;

        for (int i = 0; i < spectrum.YArray.Length; i++)
        {
            buffer[currentIndex] = spectrum.YArray[i] / window;
            double ma = 0d;
            for (int j = 0; j < window; j++)
            {
                ma += buffer[j]; 
            }

            output[i] = ma;
            currentIndex = (currentIndex + 1) % window; 
        }

        return output; 
    }
    public static IOrderedEnumerable<Tuple<int,double>> IdentifyMaxima(MzSpectrum spectrum, int movingAverageWindow, int localMaxWindowSize)
    {
        var downsampledSpectrum = DownsampleData(spectrum, 5);
        var reducedResolutionPeaks = MovingAverage(downsampledSpectrum, movingAverageWindow);
        return FindLocalMaxima(reducedResolutionPeaks, localMaxWindowSize).OrderByDescending(i => i.Item2);
    }

    public static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(List<Tuple<int,double>> indexIntensityTupleList, double[] mzValueArray, int minCharge, int maxCharge, 
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

    public static List<(double,double)> MatchChargeStateLadders(MzSpectrum spectrum, List<ChargeStateLadder> chargeStateLadders,
        List<Tuple<int, double>> indexIntensityTupleList, double ppmMatchTolerance, out List<List<int>> ladderToIndicesMap)
    {
        var mzIntensityPairs = CreateMzValueArray(spectrum, indexIntensityTupleList).OrderBy(i => i.mz).ToList(); 
        ; 
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

    public static IEnumerable<ChargeStateLadderMatch> TransformToChargeStateLadderMatch(List<List<int>> ladderToIndicesMap, List<(double,double)> mzIntensityList, List<ChargeStateLadder> ladders)
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

            List<double> chargesOfMatchingPeaks = listMzVals.Select(i => (ladder.Mass - 1.00727) / i).ToList();

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

    public static IEnumerable<Tuple<int, double>> ConvertToIndexedTuple(IEnumerable<double> array)
    {
        int index = 0; 
        foreach (double d in array)
        {
            yield return Tuple.Create(index, d);
            index++; 
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
}