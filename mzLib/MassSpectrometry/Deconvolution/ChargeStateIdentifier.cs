#nullable enable
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
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices.ComTypes;
using System.Threading;
using Accord.Collections;
using System.Threading.Tasks.Dataflow;
using Easy.Common.Interfaces;
using System.Text.RegularExpressions;
using System.Reflection;
using Easy.Common.Extensions;
using System.IO;
using System.Reflection.Metadata.Ecma335;


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
        //return DeconvoluteDataFlow(spectrumToDeconvolute, DeconvolutionParams.MinCharge, DeconvolutionParams.MaxCharge,
        //    spectrumToDeconvolute.FirstX.Value, spectrumToDeconvolute.LastX.Value,
        //    DeconvolutionParams.PeakMatchPpmTolerance, range, minimumScore: DeconvolutionParams.EnvelopeThreshold,
        //    maxThreads: DeconvolutionParams.MaxThreads).ToList(); 
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

    //public IEnumerable<IsotopicEnvelope> DeconvoluteDataFlow(MzSpectrum scan, int minCharge, int maxCharge, double minMz, double maxMz,
    //    double ppmTolerance, MzRange deconRange, int maxThreads = 19, double minimumScore = 0.6)
    //{
    //    double minimumIntensity = scan.YofPeakWithHighestY.GetValueOrDefault() * 0.01; 

    //    (int Index, double mz, double intensity)[] filteredPairs =
    //        (from i in Enumerable.Range(0, scan.XArray.Length)
    //            where scan.XArray[i] >= deconRange.Minimum && scan.XArray[i] <= deconRange.Maximum && scan.YArray[i] > minimumIntensity
    //            select (i, scan.XArray[i], scan.YArray[i]))
    //        .OrderByDescending(i => i.Item3)
    //        .ToArray();


    //    var chargeStateLadderBlock = new TransformBlock<int, IEnumerable<ChargeStateLadder>>(value =>
    //        CreateChargeStateLadders(value, scan.XArray,
    //         minCharge, maxCharge,minMz, maxMz), new ExecutionDataflowBlockOptions(){ EnsureOrdered = true, SingleProducerConstrained = true});

    //    var transformToIndexList =
    //        new TransformManyBlock<IEnumerable<ChargeStateLadder>, (ChargeStateLadder, List<int>)>(
    //            ladderList => from ladder in ladderList.AsParallel()
    //                          select (ladder, MatchChargeStateLadder(scan, ladder, ppmTolerance)), 
    //            new ExecutionDataflowBlockOptions() { MaxDegreeOfParallelism = maxThreads, EnsureOrdered = true, SingleProducerConstrained = true });

    //    var transformToLadderMatch =
    //        new TransformBlock<(ChargeStateLadder, List<int>), ChargeStateLadderMatch>(
    //            indexList =>
    //            {
    //                var match = TransformToChargeStateLadderMatch(indexList.Item2,
    //                    scan, indexList.Item1);
    //                bool success = ScoreChargeStateLadderMatch(match, scan);
                    
    //                if (!success) return null;
                    
    //                return match;
    //            }, new ExecutionDataflowBlockOptions() { MaxDegreeOfParallelism = maxThreads, EnsureOrdered = true, SingleProducerConstrained = true});

    //    var findIsotopicEnvelopesBlock = new TransformBlock<ChargeStateLadderMatch?, IEnumerable<IsotopicEnvelope>>(match =>
    //    {
    //        if (match != null)
    //        {
    //            var envelopes = FindIsotopicEnvelopes(match, scan, deconRange, minimumScore).ToList();
    //            if (envelopes.Any())
    //            {
    //                return envelopes;
    //            }
    //        }
    //        return null;
    //    }, new ExecutionDataflowBlockOptions() { MaxDegreeOfParallelism = maxThreads, EnsureOrdered = true, SingleProducerConstrained = true});


    //    var linkOptions = new DataflowLinkOptions() { PropagateCompletion = true};

    //    chargeStateLadderBlock.LinkTo(transformToIndexList, linkOptions);
    //    transformToIndexList.LinkTo(transformToLadderMatch, linkOptions); // match => match.Item1 != null && match.Item2 != null
    //    transformToLadderMatch.LinkTo(findIsotopicEnvelopesBlock, linkOptions);


    //    foreach(var pair in filteredPairs)
    //    {
    //        chargeStateLadderBlock.Post(pair.Index);
    //    }
    //    chargeStateLadderBlock.Complete();




    //    while (!findIsotopicEnvelopesBlock.Completion.IsCompleted)
    //    {
    //        bool success = findIsotopicEnvelopesBlock.TryReceive(out IEnumerable<IsotopicEnvelope> envelope); 
    //        if (success && envelope != null)
    //        {
    //            foreach (var e in envelope)
    //            {
    //                if (e.Score > minimumScore)
    //                {
    //                    yield return e;
    //                }
    //            }
    //        }
    //    }
    //}

    /// <summary>
    /// Given a spectrum and an mz range, uses the ChargeStateLadderMatch object to get the theoretical isotopic envelope, then
    /// grabs the peaks that are within the range peaks in the theoretical isotopic envelope greater than a relative intensity of 0.05. 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    /// <param name="range"></param>
    /// <returns></returns>
    internal IEnumerable<IsotopicEnvelope> FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan,
        MzRange range, double minimumScore)
    {
        ConcurrentBag<IsotopicEnvelope> results = new();
        Parallel.For(0, match.ChargesOfMatchingPeaks.Count, i =>
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
                    if (envelope.Score >= minimumScore)
                        results.Add(envelope);
                }

            }
        });
        return results;
    }

    internal IEnumerable<IsotopicEnvelope> FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan, MzRange range, 
        ConcurrentDictionary<double, IsotopicEnvelope> ieHashSet, double minimumThreshold)
    {
        this.spectrum = scan;
        Parallel.For(0, match.ChargesOfMatchingPeaks.Count, i =>
        {
            // need to be smarter than just monoisotopic mass. 
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
            
        }); 
        
        return ieHashSet.Values; 
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

    internal static IEnumerable<double> PreFilterMzVals(double[] mzVals, int minCharge, int maxCharge, 
        double minMass, double maxMass, double ppmMatchTolerance, double delta = 0.05)
    {
        double[] masses = new double[(int)((maxMass - minMass) / delta)];
        double[] counts = new double[masses.Length];

        for (int i = 0; i < masses.Length; i++)
        {
            masses[i] = minMass + delta * i;
        }

        ConcurrentDictionary<double, int> concurrentHashSet = new();
        
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
                counts[index]++;
                if (counts[index] > 1)
                {
                    concurrentHashSet.TryAdd(testMass, index);
                }
            }
        });
        
        double medianCounts = ArrayStatistics.Mean(counts.Where(i => i >= 1.1).ToArray());
        // return indexes where counts are bigger than the median. 
        HashSet<int> indexToKeepList = new();
        for (int i = 0; i < counts.Length; i++)
        {
            if (counts[i] >= medianCounts)
            {
                indexToKeepList.Add(i); 
            }
        }

        return concurrentHashSet
            .Where(i => indexToKeepList.Contains(i.Value))
            .Select(i => i.Key)
            .Distinct(new DoubleEqualityComparer());
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
    internal static IEnumerable<ChargeStateLadder> CreateChargeStateLadders(double neutralMass,
        int minCharge, int maxCharge, double minMzAllowed, double maxMzAllowed)
    {
        List<(int charge, double mz)> tempLadder = new();
        for (int j = maxCharge; j > minCharge; j--)
        {
            double mz = neutralMass.ToMz(j);
            if (mz <= maxMzAllowed && mz >= minMzAllowed)
            {
                tempLadder.Add((j, neutralMass.ToMz(j)));
            }
        }
        yield return new ChargeStateLadder(neutralMass, tempLadder.Select(k => k.mz).ToArray());
    }


    /// <summary>
    /// Produces a list of indices that match between a theoretical charge state ladder and the experimental data. 
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="chargeStateLadders"></param>
    /// <param name="ppmMatchTolerance"></param>
    /// <param name="ladderToIndicesMaps"></param>
    internal static void MatchChargeStateLaddersFast(MzSpectrum scan, List<ChargeStateLadder> chargeStateLadders,
        double ppmMatchTolerance, out List<List<int>> ladderToIndicesMaps)
    {
        var ladderDict = new ConcurrentDictionary<int, List<int>>();

        foreach(var ladder in Enumerable.Range(0, chargeStateLadders.Count))
        {
            ConcurrentBag<int> outputArray = new();
            Parallel.ForEach(Partitioner.Create(chargeStateLadders[ladder].MzVals), (t) =>
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
                        outputArray.Add(k);
                    }
                }
            });

            ladderDict.TryAdd(ladder, outputArray.ToList());
        }
        ladderToIndicesMaps = ladderDict.OrderBy(i => i.Key)
            .Select(i => i.Value)
            .ToList(); 
    }

    internal static List<int> MatchChargeStateLadder(MzSpectrum scan, ChargeStateLadder ladder, 
        double ppmMatchTolerance)
    {
        ConcurrentBag<int> output = new();
        Parallel.ForEach(ladder.MzVals, (t) =>
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
        }); 

        return output.ToList();
    }
    /// <summary>
    /// Generates a ChargeStateLadderMatch object from a charge state ladder and the list of indices that map to the original data. 
    /// </summary>
    /// <param name="ladderToIndicesMap"></param>
    /// <param name="mzIntensityList"></param>
    /// <param name="ladders"></param>
    /// <returns></returns>
    private IEnumerable<ChargeStateLadderMatch> TransformToChargeStateLadderMatch(List<List<int>> ladderToIndicesMap, 
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

            List<double> chargesOfMatchingPeaks = listMzVals.Select(i => ladder.Mass / i).ToList();

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
    /// The function that actually does the work in deconvolution, optimized for speed. 
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="deconvolutionRange"></param>
    /// <param name="spectralSimMatchThresh"></param>
    /// <returns></returns>
    private IEnumerable<IsotopicEnvelope> DeconvolutePrivateFast(MzSpectrum scan, MzRange deconvolutionRange, double spectralSimMatchThresh)
    {
        List<(double mz, double intensity)> mzIntensityPairs =
            (from i in Enumerable.Range(0, scan.XArray.Length)
             select (scan.XArray[i], scan.YArray[i])).ToList();

        // go to next if it is found in seen 

        ConcurrentBag<IsotopicEnvelope> envelopes = new();
        ConcurrentDictionary<double, IsotopicEnvelope> ieHashSet = new(new DoubleEqualityComparer());

        var neutralMasses = PreFilterMzVals(scan.XArray, DeconvolutionParams.MinCharge, DeconvolutionParams.MaxCharge, 10000, 100000,
            DeconvolutionParams.PeakMatchPpmTolerance).ToList();
        var chargeStateLadders = CreateChargeStateLadders(neutralMasses, DeconvolutionParams.MinCharge,
            DeconvolutionParams.MaxCharge, scan.FirstX.Value, scan.LastX.Value);

        Parallel.ForEach(chargeStateLadders, (ladder) =>
        {
            List<int> ladderToIndexMap = MatchChargeStateLadder(scan, ladder, DeconvolutionParams.PeakMatchPpmTolerance); 
            

            var ladderMatch = TransformToChargeStateLadderMatch(ladderToIndexMap, scan, ladder);
            var successfulMatch = ScoreChargeStateLadderMatch(ladderMatch, scan);
            if (successfulMatch)
            { 
                FindIsotopicEnvelopes(ladderMatch!, scan, deconvolutionRange, ieHashSet, DeconvolutionParams.EnvelopeThreshold);
            }
        });
        
        return ieHashSet.Values.OrderByDescending(i => i.Score); 
    }

    /// <summary>
    /// Basic function for harmonic checking, if we want to do that in the future. However, it is a little too zealous, and removes more harmonics than it should. 
    /// </summary>
    /// <param name="listEnvelopes"></param>
    /// <returns></returns>
    private IEnumerable<IsotopicEnvelope> CheckForHarmonics(IList<IsotopicEnvelope> listEnvelopes)
    {
        for (int i = 0; i < listEnvelopes.Count; i++)
        {
            bool isLowHarmonicScore = false;
            bool isLowHarmonicMass = false;
            for (int j = 0; j < listEnvelopes.Count; j++)
            {
                double scoresScore = listEnvelopes.ElementAt(i).MonoisotopicMass / listEnvelopes.ElementAt(j).MonoisotopicMass;
                isLowHarmonicScore = Math.Abs(scoresScore - 2d) <= 0.1;
                double massScore = listEnvelopes.ElementAt(i).MonoisotopicMass / listEnvelopes.ElementAt(j).MonoisotopicMass;
                isLowHarmonicMass = Math.Abs(massScore - 2d) <= 0.05;
                if (isLowHarmonicMass)
                {
                    listEnvelopes.RemoveAt(j);
                }
            }

        }
        return listEnvelopes;
    }

    public static double GetTheoreticalMostIntenseMassWithDiff(int massIndex, double diff)
    {
        return allMasses[massIndex][0] + diff; 
    }
    /// <summary>
    /// Performs spectral similarity calculation between the theoretical isotopic envelopes and the experimental data.
    /// </summary>
    /// <param name="envelope"></param>
    private void RescoreIsotopicEnvelope(IsotopicEnvelope envelope)
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

    public static double GetTheoreticalMostAbundantMass(int massIndex)
    {
        return allMasses[massIndex][0];
    }


    /// <summary>
    /// Runs scoring functions found within ChargeStateLadderMatch and performs basic filtering to remove ChargeStateLadderMatches that
    /// don't fit the data. 
    /// </summary>
    /// <param name="ladderMatches"></param>
    /// <param name="scan"></param>
    /// <returns></returns>
    private IEnumerable<ChargeStateLadderMatch?> ScoreChargeStateLadderMatches(List<ChargeStateLadderMatch> ladderMatches, MzSpectrum scan)
    {
        var orderByIntensity = ladderMatches
            //.Where(i => i.MatchingMzPeaks.Count > i.TheoreticalLadder.MzVals.Length / 3)
            .OrderByDescending(i => i.IntensitiesOfMatchingPeaks.Sum())
            .ToList();
        
        foreach (var envelope in orderByIntensity)
        {
            GetMassIndex(envelope, scan);
            List<double> monoGuesses = new();
            List<double> errorInMostIntense = new();
            //assume that the selected m / z is the most intense peak in the envelope. 

            if (envelope.TheoreticalLadder.Mass <= DeconvolutionParams.MinimumMassDa) continue;

            envelope.CompareTheoreticalNumberChargeStatesVsActual();
            if (envelope.PercentageMzValsMatched < 0.2) continue;

            envelope.CalculateChargeStateScore();
            if (envelope.SequentialChargeStateScore < -1.1) continue;

            envelope.CalculateEnvelopeScore();
            if (envelope.EnvelopeScore < 0.1) continue;

            for (int i = 0; i < envelope.MatchingMzPeaks.Count; i++)
            {
                double deChargedMass = envelope.MatchingMzPeaks[i]
                    .ToMass((int)Math.Round(envelope.ChargesOfMatchingPeaks[i]));
                // get the theoretical highest peak in isotopic envelope. 

                double monoGuess = deChargedMass - diffToMonoisotopic[envelope.MassIndex]; 

                monoGuesses.Add(monoGuess);
            }

            envelope.MonoErrors = errorInMostIntense;
            envelope.MonoGuesses = monoGuesses;
            // use median because the center of the edges of the charge state distributions are likely to be 
            // more error prone. 

            envelope.MonoisotopicMass = monoGuesses.Median();

            envelope.ScoreByIntensityExplained(scan);
            if(envelope != null) yield return envelope;
        }
    }

    private bool ScoreChargeStateLadderMatch(ChargeStateLadderMatch match, MzSpectrum scan)
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
    /// <summary>
    /// get the range of the theoretical isotopic envelopes and then pulls all the peaks within that range from the original data. 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    /// <param name="isolationRange"></param>
    /// <param name="chargeState"></param>
    /// <returns></returns>
    private IsotopicEnvelope FillIsotopicEnvelopeByBounds(ChargeStateLadderMatch match, MzSpectrum scan, MzRange isolationRange, int chargeState)
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
    // method to find a common mass within a tolerance in a list of charge state ladders, 
    // then see if the charge state ladders overlap
    internal static void ChargeStateLadderOverlap(List<ChargeStateLadder> ladderList)
    {
        // flatten the list at mass and look for duplicates
        var listMasses = ladderList.Select(i => i.Mass);
        //get counts of masses 
        var g = listMasses.GroupBy(i => i);
        foreach (var grp in g)
        {
            if (grp.Count() > 1) ; 
        }


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

public static class IsotopicEnvelopeExtensions
{
    public static IDictionary<double, IEnumerable<IsotopicEnvelope>> ObserveAdjacentChargeStates(this IEnumerable<IsotopicEnvelope> listEnvelopes)
    {
        HashSet<double> massHashSet = new();
        var orderedMonoMasses = listEnvelopes.OrderBy(i => i.MonoisotopicMass)
            .Select(i => i.MonoisotopicMass)
            .ToArray();
        var orderedArrayOfEnvelopes = listEnvelopes
            .OrderBy(i => i.MonoisotopicMass)
            .ToList();

        ConcurrentDictionary<double, IEnumerable<IsotopicEnvelope>> dictChargeStateEnvelopes = new();

        int indexer = 0;
        List<IsotopicEnvelope> tempIsotopicEnvelope = new();
        while (indexer < orderedMonoMasses.Length)
        {
            
            if (massHashSet.Contains(orderedMonoMasses[indexer]))
            {
                tempIsotopicEnvelope.Add(orderedArrayOfEnvelopes.ElementAt(indexer)); 
                indexer++; 
            }
            else
            {
                if (indexer == 0)
                {
                    massHashSet.Add(orderedMonoMasses[indexer]);
                    indexer++; 
                    continue;
                }

                dictChargeStateEnvelopes.TryAdd(orderedMonoMasses[indexer - 1], tempIsotopicEnvelope.DistinctBy(i => i.Charge));
                tempIsotopicEnvelope = new(); 

                massHashSet.Add(orderedMonoMasses[indexer]);
                tempIsotopicEnvelope.Add(orderedArrayOfEnvelopes.ElementAt(indexer));
                indexer++; 
            }

        }
        return dictChargeStateEnvelopes;
    }
}






