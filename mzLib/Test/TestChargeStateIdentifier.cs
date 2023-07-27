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
using OpenMcdf.Extensions;
using OxyPlot.Axes;


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

    [Test]
    public void TestFindIsotopicEnvelopes()
    {
        double[] mzValsShouldWork = new[] { 167.6737, 201.007, 251.007, 334.3403, 501.007, 1001.007 };
        double[] intenValsShouldWork = new[] { 0.120985, 0.176033, 0.199471, 0.176033, 0.120985, 0.064759 };
        ChargeStateLadder ladder = new(1001.007.ToMass(1), mzValsShouldWork);
        ChargeStateLadderMatch match = new ChargeStateLadderMatch();
        match.TheoreticalLadder = ladder; 
        match.IntensitiesOfMatchingPeaks = intenValsShouldWork.ToList();
        match.MatchingMzPeaks = mzValsShouldWork.ToList();
        match.ChargesOfMatchingPeaks = new List<double> { 6d, 5d, 4d, 3d, 2d, 1d };



        MzSpectrum scan = new MzSpectrum(mzValsShouldWork, intenValsShouldWork, true);
        MzRange range = new MzRange(300, 350);
        ConcurrentDictionary<double, List<IsotopicEnvelope>> isoEnvelopeHashSet = new();

        ChargeStateIdentifier csi = new(new ChargeStateDeconvolutionParams(1,2, 200, 1, 5000,20000));
        
        csi.FindIsotopicEnvelopes(match, scan, range, isoEnvelopeHashSet, 0.0);
    }

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
            { 1001.007, 953.387952380952, 910.097909090909, 870.572217391304, 834.340333333333 };

        int minCharge = 5;
        int maxCharge = 30;
        double minMass = 10000;
        double maxMass = 30000;
        double ppmMatchTolerance = 4.0;
        double delta = 1.003;
    }
}

