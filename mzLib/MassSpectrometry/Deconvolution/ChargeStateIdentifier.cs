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
using System.Runtime.InteropServices.ComTypes;
using System.Security.Cryptography.X509Certificates;
using System.Threading.Tasks;
using Accord.Statistics.Distributions.Univariate;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;
using Accord.Statistics.Testing;
using Chemistry;
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
    
    private IEnumerable<IsotopicEnvelope> DeconvoluteInternal(MzSpectrum scan, MzRange range)
    {
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
            List<ChargeStateLadderMatch?> matchList = ScoreChargeStateLadderMatches(chargeStateLadderMatches, scan).ToList();

            var bestScoringChargeStateLadderMatch = matchList
                .Where(i => i.EnvelopeScore > 0.1)
                .Where(i => i.SequentialChargeStateScore > -1.9)
                .Where(i => i.PercentageMzValsMatched > 0.2)
                .MinBy(i => Math.Abs(i.MeanError)); 

            if (bestScoringChargeStateLadderMatch == null)
            {
                iterations++;
                indexer++; 
                continue; 
            }
            //forbiddenMassIndex.Add(bestScoringChargeStateLadderMatch.MassIndex);

            var isotopicEnvelopesIntermediate = FindIsotopicEnvelopes(bestScoringChargeStateLadderMatch!, scan);

            foreach (var envelope in isotopicEnvelopesIntermediate)
            {
                if (envelope == null) continue;
                if (envelope.Peaks.Count < 2) continue;
                
                if (range == null)
                {
                    yield return envelope; 
                }
                else
                {
                    if (envelope.Peaks.Any(i => i.mz >= range.Minimum && i.mz <= range.Maximum))
                    {
                        yield return envelope;
                    }
                }
            }
            indexer++;
            iterations++;
        }
    }

    public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrumToDeconvolute, MzRange range)
    {
        return DeconvoluteInternal(spectrumToDeconvolute, range).OrderByDescending(i => i.Score);

        //return CleanUpIsotopologues(isotopicEnvelopeList, 5)
        //    .Where(i => i.Peaks.Count > 2);
        //return MergeIsotopicEnvelopes(cleanedList); 
    }

    private IEnumerable<IsotopicEnvelope> MergeIsotopicEnvelopes(IEnumerable<IsotopicEnvelope> envelopes)
    {
        var distinctMonoisotopicMasses = envelopes.Select(i => i.MonoisotopicMass).Distinct().ToArray();
        var uniqueIsotopologues = envelopes.DistinctBy(i => i.MonoisotopicMass).ToList();
        for (int i = 0; i < distinctMonoisotopicMasses.Length; i++)
        {
            uniqueIsotopologues[i].ReplacePeaksList(envelopes
                .Where(k => k.MonoisotopicMass == distinctMonoisotopicMasses[i])
                .SelectMany(k => k.Peaks)
                .Distinct()
                .ToList());
        }
        return uniqueIsotopologues;
    }

    public IEnumerable<IsotopicEnvelope> CleanUpIsotopologues(IEnumerable<IsotopicEnvelope> envelopes, double ppmTolerance)
    {
        var orderedEnvelopes = envelopes.OrderBy(i => i.MonoisotopicMass);
        foreach (var range in Enumerable.Range(0, orderedEnvelopes.Count() - 1))
        {
            if (range == 0)
            {
                yield return orderedEnvelopes.ElementAt(range);
                continue; 
            } 

            if ((orderedEnvelopes.ElementAt(range+1).MonoisotopicMass - orderedEnvelopes.ElementAt(range).MonoisotopicMass) 
                <= Constants.C13MinusC12 + ppmTolerance / 1E6 * Constants.C13MinusC12)
            {
                orderedEnvelopes.ElementAt(range + 1)
                    .SetMonoisotopicMass(orderedEnvelopes.ElementAt(range).MonoisotopicMass); 
                yield return orderedEnvelopes.ElementAt(range);
                yield return orderedEnvelopes.ElementAt(range + 1); 
            }

             
        }
    }
    public List<int> FindIsotopologues(IEnumerable<IsotopicEnvelope> envelopes, double ppmTolerance)
    {
        var orderedByMonoMass = envelopes.OrderBy(i => i.MonoisotopicMass).ToList();

        var indexList = new List<int>(); 
        for(int i = 0; i < orderedByMonoMass.Count - 1; i++)
        {
            if (Math.Abs(orderedByMonoMass[i + 1].MonoisotopicMass
                         - orderedByMonoMass[i].MonoisotopicMass) 
                <= Constants.C13MinusC12 + ppmTolerance/1E6*Constants.C13MinusC12)
            {
                indexList.Add(i);
            }
        }

        return indexList;
    }

    public IEnumerable<ChargeStateLadderMatch?> ScoreChargeStateLadderMatches(List<ChargeStateLadderMatch> ladderMatches, MzSpectrum scan)
    {

        var orderByIntensity = ladderMatches
            //.Where(i => i.MatchingMzPeaks.Count > i.TheoreticalLadder.MzVals.Length / 3)
            .OrderByDescending(i => i.IntensitiesOfMatchingPeaks.Sum())
            .ToList();
        // check for harmonics
        double[] ladderMasses = orderByIntensity.Select(i => i.TheoreticalLadder.Mass).ToArray();

        List<(int, int)> indicesList = new List<(int, int)>();
        foreach (var envelope in orderByIntensity)
        {
            GetMassIndex(envelope, scan);
            List<double> monoGuesses = new();
            List<double> errorInMostIntense = new();
            List<double> deChargedMassesList = new(); 
            // assume that the selected m/z is the most intense peak in the envelope. 
            for (int i = 0; i < envelope.MatchingMzPeaks.Count; i++)
            {
                double deChargedMass = envelope.MatchingMzPeaks[i]
                    .ToMass((int)Math.Round(envelope.ChargesOfMatchingPeaks[i]));
                deChargedMassesList.Add(deChargedMass);
                // get the theoretical highest peak in isotopic envelope. 
                double mz = allMasses[envelope.MassIndex][0]
                    .ToMz((int)Math.Round(envelope.ChargesOfMatchingPeaks[i]));
                int closestPeakMzIndex = spectrum.GetClosestPeakIndex(mz);
                double mass = spectrum.XArray[closestPeakMzIndex].ToMass((int)Math.Round(envelope.ChargesOfMatchingPeaks[i]));  
                monoGuesses.Add(mass - diffToMonoisotopic[envelope.MassIndex]);
                errorInMostIntense.Add(mass - allMasses[envelope.MassIndex][0]); 
            }
            envelope.MonoErrors = errorInMostIntense;
            envelope.MonoGuesses = monoGuesses;
            envelope.MonoisotopicMass = monoGuesses.Average(); 
            
            envelope.CalculateEnvelopeScore();
            envelope.Score = ScoreByIntensityExplained(envelope, scan); 
            envelope.PercentageMzValsMatched = CompareTheoreticalNumberChargeStatesVsActual(envelope);
            envelope.CalculateChargeStateScore();
            yield return envelope;
        }
    }

    public static double GetDiffToMonoisotopic(int massIndex)
    {
        return diffToMonoisotopic[massIndex];
    }

    //public (double rsquared, double slope) CalcRSquared(ChargeStateLadderMatch cslm)
    //{
    //    if (cslm.ChargesOfMatchingPeaks.Count < this.DeconvolutionParams.MinCharge) return (0, 0); 
    //    (double a, double b) lineParams = Fit.Line(cslm.ChargesOfMatchingPeaks.ToArray(), 
    //        Enumerable.Range(0, cslm.ChargesOfMatchingPeaks.Count)
    //            .Select(i => (double)i)
    //            .ToArray());

    //    var modeledValues = Enumerable.Range(0, cslm.ChargesOfMatchingPeaks.Count)
    //        .Select(i => (double)i * lineParams.b + lineParams.a).ToArray(); 
        
    //    return (GoodnessOfFit.CoefficientOfDetermination(modeledValues, cslm.ChargesOfMatchingPeaks), lineParams.b); 
    //}

    public double CompareTheoreticalNumberChargeStatesVsActual(ChargeStateLadderMatch ladderMatch)
    {
        return (double)ladderMatch.ChargesOfMatchingPeaks.Count / (double)ladderMatch.TheoreticalLadder.MzVals.Length;
    }
    public double ScoreByIntensityExplained(ChargeStateLadderMatch match, MzSpectrum spectrum, double threshold = 0.01)
    {
        double spectrumThresholdVal = spectrum.YArray.Max() * threshold;
        return match.IntensitiesOfMatchingPeaks.Where(i => i >= spectrumThresholdVal).Sum();
    }

    private void GetMassIndex(ChargeStateLadderMatch match, MzSpectrum scan)
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
    public double EnvelopeScore { get; set; }
    public double Score { get; set; }
    public int MassIndex { get; set; }
    public double MonoisotopicMass { get; set; }
    public double SequentialChargeStateScore { get; set; }
    public List<double> MonoGuesses { get; set; }
    public List<double> MonoErrors { get; set; }
    public double PercentageMzValsMatched { get; set; }
    public double MeanError => CalculateMonoErrorMean();

    private double CalculateMonoErrorMean()
    {
        double sum = 0; 
        for (int i = 0; i < MonoGuesses.Count; i++)
        {
            sum += MonoErrors[i]; 
        }

        return sum / MonoisotopicMass / MonoErrors.Count * 1E6; 
    }
    public void CalculateChargeStateScore()
    {
        if (!(PercentageMzValsMatched < 0.25))
        {
            var chargesList = ChargesOfMatchingPeaks
                .Zip(ChargesOfMatchingPeaks.Skip(1), (x, y) => y - x)
                .Where(i => Math.Abs(i) > 0.1);
            if (chargesList.Any())
            {
                SequentialChargeStateScore = chargesList.Average();
                return; 
            }

            SequentialChargeStateScore = -10000;
            return; 
        }
        SequentialChargeStateScore = -10000;
    }
    public (double mzMostIntense, double maxIntensityPeak) GetMostIntenseTuple()
    {

        int indexOfMostIntense = IntensitiesOfMatchingPeaks.IndexOf(IntensitiesOfMatchingPeaks.Max()); 
        return (MatchingMzPeaks[indexOfMostIntense], IntensitiesOfMatchingPeaks[indexOfMostIntense]);
    }

    internal void CalculateEnvelopeScore()
    {
        if (IntensitiesOfMatchingPeaks.Count < 4)
        {
            EnvelopeScore = 0d;
            return; 
        } 
        double[] coefficients =
            Fit.Polynomial(ChargesOfMatchingPeaks.ToArray(), IntensitiesOfMatchingPeaks.ToArray(), 2); 
        double c = coefficients[0];
        double b = coefficients[1];
        double a = coefficients[2];

        // calculate theoretical polynomial to get R^2. 
        double[] theoreticalPolynom = ChargesOfMatchingPeaks
            .Select(x => c + b*x + a*x*x)
            .ToArray();
        double sum1 = 0;
        double sum2 = 0;
        double ymean = IntensitiesOfMatchingPeaks.Mean();

        for (int i = 0; i < ChargesOfMatchingPeaks.Count; i++)
        {
            sum1 += Math.Pow(IntensitiesOfMatchingPeaks[i] - theoreticalPolynom[i], 2);
            sum2 += Math.Pow(IntensitiesOfMatchingPeaks[i] - ymean, 2); 
        }

        EnvelopeScore = 1 - sum1 / sum2;
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