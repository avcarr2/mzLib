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
using System.Windows.Input;
using Accord;
using MathNet.Numerics.Statistics;


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
        ChargeStateDeconvolutionParams deconParams = new(5, 60, 10, maxThreads: 1, envelopeThreshold:0.1);
        ChargeStateIdentifier csi = new(deconParams);

        Stopwatch watch = new();
        watch.Start();
        var results = csi.Deconvolute(scan, new MzRange(800, 2000)); 

        watch.Stop();
        Console.WriteLine(watch.ElapsedMilliseconds);
        //WriteIsotopicEnvelopesToPlottingPoints(results);
        //using var writer = new StreamWriter("deconvSpectra.csv");
        //writer.WriteLine("Mass,Intensity");
        
        //foreach (var pair in filteredIntensities)
        //{
        //    writer.WriteLine("{0},{1}", pair.Item1, pair.Item2);
        //}
        //writer.Flush();
        WriteIsotopicEnvelopesToPlottingPoints(results);
    }

    [Test]
    [TestCase(30)]
    [TestCase(100)]
    [TestCase(15)]
    [TestCase(5)]
    [TestCase(60)]
    
    public void Figure2MediumRes(double peakmatchTol)
    {
        string path = @"D:\FornelliDataSet\2015_03_19_Imr90_mh_fr_4_tech_1.raw";
        FilteringParams filteringParams = new FilteringParams(minimumAllowedIntensityRatioToBasePeak:0.1, applyTrimmingToMs1:true);
        var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList()[746];
        ChargeStateDeconvolutionParams deconParams = new(25, 60, 
            peakmatchTol, maxThreads: 1, envelopeThreshold: 0.20, 
            deltaMass:1.0, minimumMass:10000, maximumMass:150000);
        ChargeStateIdentifier csi = new(deconParams);

        Stopwatch watch = new();
        watch.Start();
        var results = csi.Deconvolute(scan.MassSpectrum, new MzRange(838, 851))
            .OrderByDescending(i => i.Score)
            .ToList(); 

        watch.Stop();

        //using var writer = new StreamWriter("deconvSpectra.csv");
        //foreach (var pair in filteredIntensities)
        //{
        //    writer.WriteLine("{0},{1}", pair.Item1, pair.Item2);
        //}
        //writer.Flush();
        Console.WriteLine(watch.ElapsedMilliseconds);
        Console.WriteLine("Total number of results: {0}", results.Count());

        WriteIsotopicEnvelopesToPlottingPoints(results, "");
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

    //[Test]
    //public void TestPreFilterMzVals()
    //{
    //    string path = @"D:\220220AveragedDatasets\LVS Jurkat\02-18-20_jurkat_td_rep2_fract9.raw";
    //    var scans = MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList();

    //    ChargeStateDeconvolutionParams deconParams = new(5, 60, 5, maxThreads: 15, envelopeThreshold: 0.6);
    //    ChargeStateIdentifier csi = new(deconParams);
    //    ConcurrentBag<int> countsOfValues = new();

    //    ConcurrentBag<IsotopicEnvelope> results = new();
    //    double[] massesArray =
    //        ChargeStateIdentifier.GenerateMassesIndex(10000, 60000, deconParams.PeakMatchPpmTolerance); 

    //    Parallel.For(0, 1, i =>
    //    {
    //        ConcurrentDictionary<double, List<IsotopicEnvelope>> ieHashSet = new(new DoubleEqualityComparer());
    //        Stopwatch watch = new();
    //        watch.Start();
    //        // slow step about 200 ms
    //        var output = ChargeStateIdentifier.PreFilterMzVals(scans[i].MassSpectrum.XArray, 
    //            scans[i].MassSpectrum.YArray, 5, 60, )
    //            .OrderBy(z => z)
    //            .ToList();
    //        // fast step
    //        var ladder = ChargeStateIdentifier.CreateChargeStateLadders(output, 5, 60, 800, 2000);

    //        foreach (var m in ladder)
    //        {
    //            var index = ChargeStateIdentifier.MatchChargeStateLadder(scans[i].MassSpectrum, m,
    //                deconParams.PeakMatchPpmTolerance);
    //            var ladderMatch =
    //                ChargeStateIdentifier.TransformToChargeStateLadderMatch(index, scans[i].MassSpectrum, m);

    //            var successfulMatch = csi.ScoreChargeStateLadderMatch(ladderMatch, scans[i].MassSpectrum);

    //            if (successfulMatch)
    //            {
    //                csi.FindIsotopicEnvelopes(ladderMatch!, scans[i].MassSpectrum, 
    //                    new MzRange(900, 1000), ieHashSet, 0.6);

    //            }
    //        }
    //        var sortedResults = ieHashSet.Values
    //            .SelectMany(i => i)
    //            .OrderByDescending(i => i.Score)
    //            .ToList(); 
    //        Console.WriteLine(watch.ElapsedMilliseconds);
    //    });
    //}

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

        return (maxPairs.mz.ToMz(envelope.Charge), maxPairs.intensity, charge, 
            envelope.MonoisotopicMass, envelope.Score, ChargeStateIdentifier.GetDiffToMonoisotopic(envelope.MassIndex)); 
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

    [Test]
    public void TestGenerateMassesIndex()
    {
        double minMass = 9500;
        double maxMass = 60000;
        double ppmTolerance = 15;

        double[] results = ChargeStateIdentifier.GenerateMassesIndex(minMass, maxMass, 15);
        Console.WriteLine("{0},{1}", results[0], results[^1]);
    }


}

