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
using MathNet.Numerics.Statistics;
using OxyPlot.Axes;
using OxyPlot.Series;
using UsefulProteomicsDatabases.Generated;


namespace Test;
public class TestChargeStateIdentifier
{
    //[Test]
    //public void TestCreateChargeStateLadders()
    //{
    //    ChargeStateDeconvolutionParams deconParams = new(1, 5, 5, maxThreads:12);
    //    ChargeStateIdentifier csi = new(deconParams);
    //    double[] testMzVals = new[] { 1000.0 }; 

    //    IEnumerable<ChargeStateLadder> ladders = ChargeStateIdentifier.CreateChargeStateLadders(0, testMzVals, 
    //        deconParams.MinCharge, deconParams.MaxCharge, 
    //        500d, 2000d);
    //    var massArrayTest = ladders.Select(i => i.Mass).ToArray();
    //    var massArrayExpected = new[] { 998.9927, 1997.98545, 2996.9782, 3995.97089, 4994.96361 };
    //    Assert.That(massArrayExpected, Is.EqualTo(massArrayTest).Within(0.01));
    //    var valsArray = ladders.Select(i => i.MzVals.Length).ToArray();
    //    var expectedArrayLengths = new int[] { 2, 4, 4, 4, 3 };
    //    Assert.That(valsArray, Is.EqualTo(expectedArrayLengths));
    //}

    //[Test]
    //public void TestChargeStateLadderMatch()
    //{
    //    double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007};
    //    double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759}; 

    //    MzSpectrum spectrumShouldWork = new MzSpectrum(mzValsShouldWork, intenValsShouldWork, true);

    //    ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork); 
    //    ChargeStateLadderMatch match = new();
    //    match.TheoreticalLadder = ladder; 
    //    match.IntensitiesOfMatchingPeaks = intenValsShouldWork.ToList(); 
    //    match.MatchingMzPeaks = mzValsShouldWork.ToList();
    //    match.ChargesOfMatchingPeaks = Enumerable.Range(1, 6).Reverse().Select(i =>(double)i).ToList();

    //    match.ScoreByIntensityExplained(spectrumShouldWork, threshold:0); 
    //    match.CalculateEnvelopeScore();
    //    match.CalculateChargeStateScore();
    //    match.CompareTheoreticalNumberChargeStatesVsActual(); 

    //    Assert.That(match.EnvelopeScore, Is.EqualTo(0.9765).Within(0.05));
    //    Assert.That(match.Score, Is.EqualTo(0.858).Within(0.05));
    //    Assert.That(match.PercentageMzValsMatched, Is.EqualTo(1.0).Within(0.01));
    //    Assert.That(match.SequentialChargeStateScore, Is.EqualTo(-1d));
    //}

    //[Test]
    //public void TestChargeStateLadderFailures()
    //{
    //    double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007 };
    //    double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759 };
    //    ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork);

    //    double[] mzValsShouldntWork = new[] { 201.007, 334.3403, 1001.007 };
    //    double[] intenValsShouldntWork = new[] { 0.176033, 0.176033, 0.064759 };

    //    MzSpectrum spectrumShouldWork = new MzSpectrum(mzValsShouldntWork, intenValsShouldntWork, true);

    //    ChargeStateLadderMatch match = new();
    //    match.TheoreticalLadder = ladder;
    //    match.IntensitiesOfMatchingPeaks = intenValsShouldntWork.ToList();
    //    match.MatchingMzPeaks = mzValsShouldntWork.ToList();
    //    match.ChargesOfMatchingPeaks = new List<double>() { 6, 4, 2 }; 
        
    //    match.ScoreByIntensityExplained(spectrumShouldWork, threshold: 0);
    //    match.CalculateEnvelopeScore();
    //    match.CalculateChargeStateScore();
    //    match.CompareTheoreticalNumberChargeStatesVsActual();

    //    Assert.That(match.EnvelopeScore, Is.EqualTo(0d));
    //    Assert.That(match.Score, Is.EqualTo(intenValsShouldntWork.Sum()));
    //    Assert.That(match.PercentageMzValsMatched, Is.EqualTo(0.5).Within(0.01));
    //    Assert.That(match.SequentialChargeStateScore, Is.EqualTo(-2d));
    //}

    //[Test]
    //public void TestFindIsotopicEnvelopes()
    //{
    //    double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007 };
    //    double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759 };
    //    ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork);
    //    ChargeStateLadderMatch match = new ChargeStateLadderMatch();
    //    match.TheoreticalLadder = ladder; 
    //    match.IntensitiesOfMatchingPeaks = intenValsShouldWork.ToList();
    //    match.MatchingMzPeaks = mzValsShouldWork.ToList();
    //    match.ChargesOfMatchingPeaks = new List<double> { 6d, 5d, 4d, 3d, 2d, 1d };



    //    MzSpectrum scan = new MzSpectrum(mzValsShouldWork, intenValsShouldWork, true);
    //    MzRange range = new MzRange(300, 350);
    //    ConcurrentDictionary<double, List<IsotopicEnvelope>> isoEnvelopeHashSet = new();

    //    ChargeStateIdentifier csi = new(new ChargeStateDeconvolutionParams(1,2, 200, 1, 5000,20000));
        
    //    csi.FindIsotopicEnvelopes(match, scan, range, isoEnvelopeHashSet, 0.0);
    //}

    //[Test]
    //[TestCase(2,60)]
    //[Repeat(5)]
    //public void TestDeconvolution(int minCharge, int maxCharge)
    //{
    //    string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
    //        "chargeStateDeconTest3raw.raw");
    //    FilteringParams filteringParams = new FilteringParams();
    //    var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList().First().MassSpectrum;

    //    ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 5, maxThreads:19); 
    //    ChargeStateIdentifier csi = new(deconParams);

    //    Stopwatch watch = new(); 
    //    watch.Start();
    //    var results = csi.Deconvolute(scan, new MzRange(600, 640)).ToList();
    //    watch.Stop();
    //    Console.WriteLine(watch.ElapsedMilliseconds);
    //}
    
    //[Test]
    //public void TestNonUniqueChargeStateLadders()
    //{
    //    double[] testMzAxis = new[]
    //        { 1001.007,953.387952380952,910.097909090909,870.572217391304,834.340333333333 };
    //    double[] testIntAxis = new[] {
    //        1000d, 2500d, 3500d, 2500d, 1000d };
    //    MzSpectrum testSpectra = new MzSpectrum(testMzAxis, testIntAxis, true);
    //    ChargeStateDeconvolutionParams deconParams = new(5, 30, 5, maxThreads: 15);
    //    ChargeStateIdentifier csi = new(deconParams);

    //    List<List<ChargeStateLadder>> laddersList = new(); 
    //    for (int i = 0; i < testMzAxis.Length; i++)
    //    {
    //        var tempList = ChargeStateIdentifier.CreateChargeStateLadders(i, testMzAxis, deconParams.MinCharge, deconParams.MaxCharge,
    //            500, 2000).ToList(); 
    //        laddersList.Add(tempList);
    //    }

    //    var masses = laddersList.Select(i => i.Select(j => j.Mass));
    //    foreach (var mass in masses.GroupBy(i => i))
    //    {
    //        if (mass.Count() > 1)
    //        {
    //            Console.WriteLine(mass);
    //        }
    //    }
    //}

    //[Test]
    //public void TestPreFiltering()
    //{
    //    double[] testMzAxis = new[]
    //        { 1001.007, 953.387952380952, 910.097909090909, 870.572217391304, 834.340333333333 };

    //    int minCharge = 5;
    //    int maxCharge = 30;
    //    double minMass = 10000;
    //    double maxMass = 30000;
    //    double ppmMatchTolerance = 4.0;
    //    double delta = 1.003;
    //}

    //[Test]
    //public void TestGlyceraldehyde()
    //{
    //    string path = @"C:\Xcalibur\data\glyceraldehyde_dehydrogenase_fornelliDataset.raw";

    //    var spectrum = Readers.ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList().First();  

    //    ChargeStateIdentifier csi = new ChargeStateIdentifier(
    //        new ChargeStateDeconvolutionParams(5, 60, 100, 19, 
    //            envelopeThreshold: 0.0, sequentialChargeStateDiff: 0.8, envelopeScoreThresh:0.6, percentageMatchedThresh:0d,
    //            deltaMass:0.1, deconType: PreFilterDeconvolutionType.Multiplicative));

    //    var result = csi.Deconvolute(spectrum.MassSpectrum, new MzRange(500, 2000)).ToList(); 

    //}

    [Test]
    public void TestIonTrap()
    {
        string path = @"D:\DeconvolutionPaper\SixProtStandardMixMeth1-firstPeak.raw"; 

        var spectrum = Readers.ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList().First();

        ChargeStateIdentifier csi = new ChargeStateIdentifier(
            new ChargeStateDeconvolutionParams(5, 60,
                10, 19, 9000, maximumMass: 60000, 
                deltaMass:0.1, envelopeThreshold: 0.2, sequentialChargeStateDiff: 0.9, 
                envelopeScoreThresh: 0.6, percentageMatchedThresh: 0.5, 
                deconType: PreFilterDeconvolutionType.Multiplicative));

        var result = csi.Deconvolute(spectrum.MassSpectrum, new MzRange(800, 860)).ToList();
        
    }

    [Test]
    public void TestWhyHarmonic()
    {
        string path = @"D:\DeconvolutionPaper\SixProtStandardMixMeth1.mzML";
        var spectrum = Readers.MsDataFileReader.GetDataFile(path);
        FilteringParams filteringParams = new FilteringParams(minimumAllowedIntensityRatioToBasePeak:0.01, applyTrimmingToMs1:true);
        spectrum.InitiateDynamicConnection();
        var scan = spectrum.GetOneBasedScanFromDynamicConnection(1469, filteringParams); 

        ChargeStateIdentifier csi = new ChargeStateIdentifier(
            new ChargeStateDeconvolutionParams(5, 100,
                15, 19, 9000, maximumMass: 70000,
                deltaMass: 0.25, envelopeThreshold: 0.01, sequentialChargeStateDiff: 0.9,
                envelopeScoreThresh: 0.6, percentageMatchedThresh: 0.5,
                deconType: PreFilterDeconvolutionType.Multiplicative));

        var results = csi.Deconvolute(scan.MassSpectrum, 
            new MzRange(scan.MassSpectrum.FirstX.Value, scan.MassSpectrum.LastX.Value))
            .DistinctBy(i => i.MonoisotopicMass)
            .ToList(); 
        spectrum.CloseDynamicConnection();

    }

    [Test]
    [TestCase(@"D:\AA_08-01-23_ISF\SixProtStandardMixMeth1.raw")]
    public void DeconProteinAGAndCytochromeC(string path)
    {
        string proteinAgPath = @"D:\AA_08-01-23_ISF\proteinAg.txt";
        string cytoCPath = @"D:\AA_08-01-23_ISF\cytoC.txt"; 

        int proteinAgScanNumber = 2244;
        int cytocScanNumber = 2523;

        using var streamWriterAg = new StreamWriter(proteinAgPath, false);
        using var cytoCyWriter = new StreamWriter(cytoCPath, false);

        var spectrum = MsDataFileReader.GetDataFile(path);
        FilteringParams filteringParams = new FilteringParams(minimumAllowedIntensityRatioToBasePeak: 0.1, applyTrimmingToMs1: true);
        spectrum.InitiateDynamicConnection();
        var scanAG = spectrum.GetOneBasedScanFromDynamicConnection(2244, filteringParams);
        var cytoScan = spectrum.GetOneBasedScanFromDynamicConnection(2523, filteringParams); 

        ChargeStateIdentifier csi = new ChargeStateIdentifier(
            new ChargeStateDeconvolutionParams(5, 100,
                15, 19, 9000, maximumMass: 70000,
                deltaMass: 0.01, envelopeThreshold: 0.01, sequentialChargeStateDiff: 0.9,
                envelopeScoreThresh: 0.6, percentageMatchedThresh: 0.5,
                deconType: PreFilterDeconvolutionType.Multiplicative));

        var resultsAg = csi.Deconvolute(scanAG.MassSpectrum,
                new MzRange(scanAG.MassSpectrum.FirstX.Value, scanAG.MassSpectrum.LastX.Value))
            .ToList();
        var resultsCyto = csi.Deconvolute(cytoScan.MassSpectrum,
            new MzRange(cytoScan.MassSpectrum.FirstX.Value, cytoScan.MassSpectrum.LastX.Value))
            .ToList();

        foreach (var i in resultsAg)
        {
            streamWriterAg.WriteLine(string.Join("\t", i.Charge, i.MonoisotopicMass));
        }
        streamWriterAg.Flush();

        foreach (var i in resultsCyto)
        {
            cytoCyWriter.WriteLine(string.Join("\t", i.Charge, i.MonoisotopicMass));
        }
        cytoCyWriter.Flush();
        spectrum.CloseDynamicConnection();
    }

}

