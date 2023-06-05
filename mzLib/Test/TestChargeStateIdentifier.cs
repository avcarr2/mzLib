using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using Chemistry;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MzLibUtil;
using Readers;
using NUnit.Framework; 
using MassSpectrometry;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Documents;
using Accord;


namespace Test;
public class TestChargeStateIdentifier
{
    [Test]
    public void TestCreateChargeStateLadders()
    {
        ChargeStateDeconvolutionParams deconParams = new(1, 5, 5, maxThreads:12);
        ChargeStateIdentifier csi = new(deconParams);
        double[] testMzVals = new[] { 1000.0 }; 

        IEnumerable<ChargeStateLadder> ladders = ChargeStateIdentifier.CreateChargeStateLadders(0, testMzVals, 
            deconParams.MinCharge, deconParams.MaxCharge, 
            500d, 2000d);
        var massArrayTest = ladders.Select(i => i.Mass).ToArray();
        var massArrayExpected = new[] { 998.9927, 1997.98545, 2996.9782, 3995.97089, 4994.96361 };
        Assert.That(massArrayExpected, Is.EqualTo(massArrayTest).Within(0.01));
        var valsArray = ladders.Select(i => i.MzVals.Length).ToArray();
        var expectedArrayLengths = new int[] { 2, 4, 4, 4, 3 };
        Assert.That(valsArray, Is.EqualTo(expectedArrayLengths));
    }

    [Test]
    public void TestChargeStateLadderMatch()
    {
        double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007};
        double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759}; 

        MzSpectrum spectrumShouldWork = new MzSpectrum(mzValsShouldWork, intenValsShouldWork, true);

        ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork); 
        ChargeStateLadderMatch match = new();
        match.TheoreticalLadder = ladder; 
        match.IntensitiesOfMatchingPeaks = intenValsShouldWork.ToList(); 
        match.MatchingMzPeaks = mzValsShouldWork.ToList();
        match.ChargesOfMatchingPeaks = Enumerable.Range(1, 6).Reverse().Select(i =>(double)i).ToList();

        match.ScoreByIntensityExplained(spectrumShouldWork, threshold:0); 
        match.CalculateEnvelopeScore();
        match.CalculateChargeStateScore();
        match.CompareTheoreticalNumberChargeStatesVsActual(); 

        Assert.That(match.EnvelopeScore, Is.EqualTo(0.9765).Within(0.05));
        Assert.That(match.Score, Is.EqualTo(0.858).Within(0.05));
        Assert.That(match.PercentageMzValsMatched, Is.EqualTo(1.0).Within(0.01));
        Assert.That(match.SequentialChargeStateScore, Is.EqualTo(-1d));
    }

    [Test]
    public void TestChargeStateLadderFailures()
    {
        double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007 };
        double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759 };
        ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork);

        double[] mzValsShouldntWork = new[] { 201.007, 334.3403, 1001.007 };
        double[] intenValsShouldntWork = new[] { 0.176033, 0.176033, 0.064759 };

        MzSpectrum spectrumShouldWork = new MzSpectrum(mzValsShouldntWork, intenValsShouldntWork, true);

        ChargeStateLadderMatch match = new();
        match.TheoreticalLadder = ladder;
        match.IntensitiesOfMatchingPeaks = intenValsShouldntWork.ToList();
        match.MatchingMzPeaks = mzValsShouldntWork.ToList();
        match.ChargesOfMatchingPeaks = new List<double>() { 6, 4, 2 }; 
        
        match.ScoreByIntensityExplained(spectrumShouldWork, threshold: 0);
        match.CalculateEnvelopeScore();
        match.CalculateChargeStateScore();
        match.CompareTheoreticalNumberChargeStatesVsActual();

        Assert.That(match.EnvelopeScore, Is.EqualTo(0d));
        Assert.That(match.Score, Is.EqualTo(intenValsShouldntWork.Sum()));
        Assert.That(match.PercentageMzValsMatched, Is.EqualTo(0.5).Within(0.01));
        Assert.That(match.SequentialChargeStateScore, Is.EqualTo(-2d));
    }
    
    // testing the scoring function is going to be a pain in the ass. 

    [Test]
    [TestCase(2,60)]
    [Repeat(5)]
    public void TestDeconvolution(int minCharge, int maxCharge)
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "chargeStateDeconTest3raw.raw");
        FilteringParams filteringParams = new FilteringParams();
        var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList().First().MassSpectrum;

        ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 5, maxThreads:19); 
        ChargeStateIdentifier csi = new(deconParams);

        Stopwatch watch = new(); 
        watch.Start();
        var results = csi.Deconvolute(scan, new MzRange(600, 640)).ToList();
        watch.Stop();
        Console.WriteLine(watch.ElapsedMilliseconds);
    }

    //[Test]
    //public void TestGroundTruthDataSet()
    //{
    //    string path = Path.Combine(@"C:\Users\Austin\Documents\Projects\MsDataSimulatorOutput\train.mzML");
    //    MsDataFile file = MsDataFileReader.GetDataFile(path); 
    //    file.InitiateDynamicConnection();
    //    MsDataScan scan = file.GetOneBasedScanFromDynamicConnection(116);
    //    ChargeStateDeconvolutionParams deconParams = new ChargeStateDeconvolutionParams(5, 120, 1, maxThreads: 12);
    //    ChargeStateIdentifier decon = new ChargeStateIdentifier(deconParams);
    //    var results = decon.Deconvolute(scan.MassSpectrum, new MzRange(900d, 1000d))
    //        .ToList();
    //    var chargeStateEnvelopes = results
    //        .ObserveAdjacentChargeStates()
    //        .OrderByDescending(i => i.Key);
        
    //    // you could score by euclidean distance and then cluster the peaks together i think. 
    //}

    // [Test]
    // public void TestDeconvolutionBigThings()
    // {
    //     
    //     int minCharge = 5;
    //     int maxCharge = 100; 
    //     string path = @"D:\DeconDataSet\SEC4-08AUG16_5uLinj_3SEC_000021 (2)..mzML";
    //     var reader = MsDataFileReader.GetDataFile(path); 
    //     reader.InitiateDynamicConnection();
    //     var scan = reader.GetOneBasedScanFromDynamicConnection(679).MassSpectrum;
    //
    //     ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 20);
    //     ChargeStateIdentifier csi = new(deconParams);
    //
    //     var results = csi.Deconvolute(scan, new MzRange(600,650)).ToList();
    //     
    //     results
    //        // .CleanUpEnvelopes()
    //         .OrderByDescending(i => i.Score)
    //         //.DistinctBy(i => i.MonoisotopicMass)
    //         .ForEach(i =>
    //         {
    //             if (i != null)
    //             {
    //                 Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
    //             }
    //         });
    //
    //     // get the charge states of the isotopic envelopes
    //     var dictChargeStateEnvelopes = GetEnvelopes(results)
    //         // possible criteria for removing: 
    //         // non-consecutive charge states; sum of peaks explain less than a certain threshold; 
    //         // peaks have less
    //         .OrderByDescending(i => i.Value.Sum(j => j.Score))
    //         .Take(30)
    //         .ToDictionary(i => i.Key, i => i.Value);
    //     // to plot: 
    // }

    // [Test]
    // public void VeryBigProteinDeconvolutionTest()
    // {
    //     int minCharge = 50;
    //     int maxCharge = 250;
    //     string path = @"D:\DeconDataSet\SEC4-08AUG16_5uLinj_3SEC_000021 (2)..mzML";
    //     var reader = MsDataFileReader.GetDataFile(path);
    //     reader.InitiateDynamicConnection();
    //     var scan = reader.GetOneBasedScanFromDynamicConnection(1016).MassSpectrum;
    //
    //     ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 10);
    //     ChargeStateIdentifier csi = new(deconParams);
    //
    //     var results = csi.Deconvolute(scan, new MzRange(985, 1000)).ToList();
    //
    //     results
    //         // .CleanUpEnvelopes()
    //         .OrderByDescending(i => i.TotalIntensity)
    //         //.DistinctBy(i => i.MonoisotopicMass)
    //         .ForEach(i =>
    //         {
    //             if (i != null)
    //             {
    //                 Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
    //             }
    //         });
    //     // suppose that you have isolated a 
    //
    // }

    //private MzSpectrum AddMzSpectraTogether(IEnumerable<MzSpectrum> spectra, double tolerancePpm)
    //{
    //    Dictionary<double, double> spectraDict = new(); 
    //    List<double> mzValues = new List<double>();
    //    List<double> intensityValue = new List<double>();

    //    foreach (var spectrum in spectra)
    //    {
    //        for (int i = 0; i < spectrum.XArray.Length; i++)
    //        {
    //            if (spectraDict.ContainsKey(spectrum.XArray[i]))
    //            {
    //                spectraDict[spectrum.XArray[i]] += spectrum.YArray[i]; 
    //            }
    //            else
    //            {
    //                spectraDict.Add(spectrum.XArray[i], spectrum.YArray[i]);
    //            }
    //        }
    //    }



    //    return new MzSpectrum(spectraDict.Select(i => i.Key)
    //            .ToArray(),
    //        spectraDict.Select(i => i.Value).ToArray(), true); 
    //}


    public List<(double xVal, double yVal, int chargeState, double monoisotopicMass)> GetMzChargeStateIntensityMonoisotopicMass(Dictionary<double, List<IsotopicEnvelope>> dictEnvelope)
    {
        List<(double xVal, double yVal, int chargeState, double monoisotopicMass)> outputList =
            new List<(double xVal, double yVal, int chargeState, double monoisotopicMass)>();

        foreach (var kvp in dictEnvelope)
        {
            double mass = kvp.Key;
            
            foreach (var envelope in kvp.Value)
            {
                double mostAbundantMz =
                    (ChargeStateIdentifier.GetDiffToMonoisotopic(envelope.MassIndex) + mass).ToMz(envelope.Charge);
                // retrieve the intensity corresponding to the mostabundant isotopic mass
                var arrayOfMzVals = envelope.Peaks
                    .OrderBy(i => i.mz)
                    .ToArray(); 

                int indexOfMostAbundant = Array.BinarySearch(arrayOfMzVals.Select(i => i.mz).ToArray(), mostAbundantMz); 

                indexOfMostAbundant = indexOfMostAbundant < 0 ? ~indexOfMostAbundant : indexOfMostAbundant;
                indexOfMostAbundant = indexOfMostAbundant >= envelope.Peaks.Count ? envelope.Peaks.Count: indexOfMostAbundant;
                var maxPairs = arrayOfMzVals[indexOfMostAbundant];
                

                int charge = envelope.Charge;
                outputList.Add((maxPairs.mz, maxPairs.intensity, charge, mass));
            }
        }
        return outputList;
    }

    [Test]
    [Repeat(1)]
    public void Figure1TestCaseEasilyResolved()
    {
        string path = @"D:\MSV000084001_StandardProteins\190226_FIlg_3_FD_500ng-averaged.raw";
        FilteringParams filteringParams = new FilteringParams();
        var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList().First().MassSpectrum;
        ChargeStateDeconvolutionParams deconParams = new(5, 60, 5, maxThreads: 1, envelopeThreshold:0.6);
        ChargeStateIdentifier csi = new(deconParams);

        Stopwatch watch = new();
        watch.Start();
        var results = csi.Deconvolute(scan, new MzRange(895,905)).OrderByDescending(i => i.Score).ToList();
        watch.Stop();
        Console.WriteLine(watch.ElapsedMilliseconds);
        //WriteIsotopicEnvelopesToPlottingPoints(results);
    }

    [Test]
    public void TestTimingWholeFile()
    {
        string path = @"D:\220220AveragedDatasets\LVS Jurkat\02-18-20_jurkat_td_rep2_fract9.raw";
        FilteringParams filteringParams = new FilteringParams(minimumAllowedIntensityRatioToBasePeak:0.1, applyTrimmingToMs1:true, applyTrimmingToMsMs:false);

        var scans = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList().Skip(1200).Take(100).ToList();
        ChargeStateDeconvolutionParams deconParams = new(5, 60, 5, maxThreads: 15, envelopeThreshold: 0.6);
        ChargeStateIdentifier csi = new(deconParams);
        ConcurrentBag<double> times = new(); 
        Parallel.ForEach(scans, (scan) =>
        {
            if (scan.MsnOrder == 1)
            {
                Stopwatch watch = new();
                watch.Start();
                var deconResults = csi.Deconvolute(scan.MassSpectrum, new MzRange(900, 1008)).ToList();
                watch.Stop();
                times.Add(watch.ElapsedMilliseconds);
            }
        });
        times.ToList().ForEach(Console.WriteLine);
    }

    [Test]
    public void DoubleEqualityComparerTest()
    {
        var double1 = 18793.58261;
        var double2 = 18793.58329; 

        var test = Math.Round(double1,2) == Math.Round(double2,2);
        Assert.True(test);

        DoubleEqualityComparer comparer = new();
        Assert.True(comparer.Equals(double1, double2)); 
    }

    [Test]
    public void TestPreFilterMzVals()
    {
        string path = @"D:\220220AveragedDatasets\LVS Jurkat\02-18-20_jurkat_td_rep2_fract9.raw";
        var scans = MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList();

        ChargeStateDeconvolutionParams deconParams = new(5, 60, 5, maxThreads: 15, envelopeThreshold: 0.6);
        ChargeStateIdentifier csi = new(deconParams);
        ConcurrentBag<int> countsOfValues = new();

        ConcurrentBag<IsotopicEnvelope> results = new(); 
        Parallel.For(0, 1, i =>
        {
            ConcurrentDictionary<double, IsotopicEnvelope> ieHashSet = new(new DoubleEqualityComparer());
            Stopwatch watch = new();
            watch.Start();
            // slow step about 200 ms
            var output = ChargeStateIdentifier.PreFilterMzVals(scans[i].MassSpectrum.XArray, 
                scans[i].MassSpectrum.YArray, 5, 60, 10000, 
                60000, 5d, 1.003)
                .OrderBy(z => z)
                .ToList();
            // fast step
            var ladder = ChargeStateIdentifier.CreateChargeStateLadders(output, 5, 60, 800, 2000);

            foreach (var m in ladder)
            {
                var index = ChargeStateIdentifier.MatchChargeStateLadder(scans[i].MassSpectrum, m,
                    deconParams.PeakMatchPpmTolerance);
                var ladderMatch =
                    ChargeStateIdentifier.TransformToChargeStateLadderMatch(index, scans[i].MassSpectrum, m);

                var successfulMatch = csi.ScoreChargeStateLadderMatch(ladderMatch, scans[i].MassSpectrum);

                if (successfulMatch)
                {
                    csi.FindIsotopicEnvelopes(ladderMatch!, scans[i].MassSpectrum, 
                        new MzRange(900, 1000), ieHashSet, 0.6);

                }
            }
            var sortedResults = ieHashSet.Values.OrderByDescending(i => i.Score).ToList(); 
            Console.WriteLine(watch.ElapsedMilliseconds);
        });
    }

    public void WriteIsotopicEnvelopesToPlottingPoints(IEnumerable<IsotopicEnvelope> envelopeCollection, string extension = ".tsv")
    {
        if (extension == ".tsv")
        {
            using var writer = new StreamWriter("plottingPoints.tsv");
            writer.WriteLine("mz\tintensity\tchargeState\tmonoMass\tscore\tdiff");
            foreach (var envelope in envelopeCollection)
            {
                var info = GetPlottingTuple(envelope); 
                writer.WriteLine(string.Join("\t", info.xVal, info.yVal, info.chargeState, info.monoMass, info.score, info.diff));
            }

            writer.Flush(); 
        }

        if (extension == "")
        {
            Console.WriteLine("mz\tintensity\tchargeState\tmonoMass\tscore\tdiff");
            foreach (var envelope in envelopeCollection)
            {
                Console.WriteLine(string.Join("\t", GetPlottingTuple(envelope)));
            }
        }
    }

    public (double xVal, double yVal, int chargeState, double monoMass, double score, double diff) GetPlottingTuple(IsotopicEnvelope envelope)
    {
        double mostAbundantMz = envelope.MonoisotopicMass + ChargeStateIdentifier.GetDiffToMonoisotopic(envelope.MassIndex); 
        // retrieve the intensity corresponding to the mostabundant isotopic mass
        double[] arrayOfMzVals = envelope.Peaks
            .OrderBy(i => i.mz)
            .Select(i => i.mz)
            .ToArray();

        int indexOfMostAbundant = Array.BinarySearch(arrayOfMzVals, mostAbundantMz);

        int indexOfMostAbundant2 = indexOfMostAbundant < 0 ? ~indexOfMostAbundant : indexOfMostAbundant;
        int indexOfMostAbundant3 = indexOfMostAbundant2 >= envelope.Peaks.Count ? envelope.Peaks.Count - 1 : indexOfMostAbundant2;
        // need to compare to see which index the x value is closer to. false if need to decrement. 
        if (indexOfMostAbundant3 > 0)
        {
            bool roundUp = arrayOfMzVals[indexOfMostAbundant3] - mostAbundantMz < mostAbundantMz - arrayOfMzVals[indexOfMostAbundant3 - 1];
            indexOfMostAbundant3 = roundUp ? indexOfMostAbundant3 : indexOfMostAbundant3 - 1;
        }
        
        // the OrderBy shouldn't be repeated, but I'm tired and about to go home. 
        var maxPairs = envelope.Peaks.OrderBy(i => i.mz).ToArray()[indexOfMostAbundant3];

        int charge = envelope.Charge;
        double mz = maxPairs.mz.ToMz(envelope.Charge); 

        return (maxPairs.mz.ToMz(envelope.Charge), maxPairs.intensity, charge, envelope.MonoisotopicMass, envelope.Score, ChargeStateIdentifier.GetDiffToMonoisotopic(envelope.MassIndex)); 
    }

    public Dictionary<double,List<IsotopicEnvelope>> GetEnvelopes(IEnumerable<IsotopicEnvelope> envelopes)
    {
        // a charge state envelope shares a monoisotopic mass 
        Dictionary<double, List<IsotopicEnvelope>> resultsDict = new(); 
        // get the list of distinct monoisotopic masses. 
        var distinctMonoisotopicMasses = envelopes
            .Select(i => i.MonoisotopicMass)
            .Distinct();

        foreach (var monoMass in distinctMonoisotopicMasses)
        {
            // retrieve all the same monoisotopic masses 
            var envelopesList = envelopes
                .Where(i => i.MonoisotopicMass == monoMass).ToList(); 

            resultsDict.Add(monoMass, envelopesList);
        }
        return resultsDict;
    }

    [Test]
    public void TestNonUniqueChargeStateLadders()
    {
        double[] testMzAxis = new[]
            { 1001.007,953.387952380952,910.097909090909,870.572217391304,834.340333333333 };
        double[] testIntAxis = new[] {
            1000d, 2500d, 3500d, 2500d, 1000d };
        MzSpectrum testSpectra = new MzSpectrum(testMzAxis, testIntAxis, true);
        ChargeStateDeconvolutionParams deconParams = new(5, 30, 5, maxThreads: 15);
        ChargeStateIdentifier csi = new(deconParams);

        List<List<ChargeStateLadder>> laddersList = new(); 
        for (int i = 0; i < testMzAxis.Length; i++)
        {
            var tempList = ChargeStateIdentifier.CreateChargeStateLadders(i, testMzAxis, deconParams.MinCharge, deconParams.MaxCharge,
                500, 2000).ToList(); 
            laddersList.Add(tempList);
        }

        var masses = laddersList.Select(i => i.Select(j => j.Mass));
        foreach (var mass in masses.GroupBy(i => i))
        {
            if (mass.Count() > 1)
            {
                Console.WriteLine(mass);
            }
        }
    }

    [Test]
    public void TestPreFiltering()
    {
        double[] testMzAxis = new[]
            { 1001.007,953.387952380952,910.097909090909,870.572217391304,834.340333333333 };

        int minCharge = 5;
        int maxCharge = 30;
        double minMass = 10000;
        double maxMass = 30000; 
        double ppmMatchTolerance = 4.0;
        double delta = 0.5; 

        List<(int, double)> massesList = new();

        double[] masses = new double[(int)((maxMass - minMass) / delta)];
        int[] counts = new int[masses.Length];

        for (int i = 0; i < masses.Length; i++)
        {
            masses[i] = minMass + delta * i;
        }

        List<double> neutralMasses = new(); 

        for (int i = maxCharge; i >= minCharge; i--)
        {
            for (int j = 0; j < testMzAxis.Length; j++)
            {
                var testMass = testMzAxis[j].ToMass(i);
                if (testMass > maxMass || testMass < minMass)
                {
                    continue; 
                }

                int index = GetBucket(masses, testMass);
                counts[index]++;

                if (counts[index] > 1)
                {
                    if (!neutralMasses.Any(d => Math.Abs(testMass - d) * 1e6 / testMass <= ppmMatchTolerance))
                    {
                        neutralMasses.Add(testMass);
                    }
                }
            }
        }
        


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

    private static bool DoubleAreNearlyEqual(double d1, double d2, double ppmMatchTolerance)
    {
        return Math.Abs(d1 - d2) / 1e6 * d1 < ppmMatchTolerance;
    }


}

