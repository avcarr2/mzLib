using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using Chemistry;
using Easy.Common.Extensions;
using MzLibUtil;
using Readers;
using NUnit.Framework; 
using MassSpectrometry;


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

        //var chargeState = results.ToList()
        //    .ObserveAdjacentChargeStates()
        //    .OrderByDescending(i => i.Key)
        //    //.DistinctBy(i => i.MonoisotopicMass)
        //    ;


        //char[] separators = new char[] { '(', ')', ','}; 
        //Console.WriteLine("mz\tintensity\tchargeState\tmonoMass\tscore");
        //List<(double, double, int, double, double, double)> plottingTuples = new List<(double, double, int, double, double, double)>();

        //chargeState
        //    .Where(j => j.Value.Count() > 2)
        //    .ForEach(j =>
        //{
        //    foreach (var envelope in j.Value)
        //    {
        //        plottingTuples.Add(GetPlottingTuple(envelope));
        //        //Console.WriteLine(string.Join("\t",GetPlottingTuple(envelope)
        //        //    .ToString()
        //        //    .Split(separators))); 
        //    }

        //    //Console.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}",
        //    //    j.Key,
        //    //    string.Join(";", j.Value.Select(i => i.Charge)),
        //    //    string.Join(",", string.Join(",", j.Value.Select(i => i.Score))),
        //    //    string.Join(",", j.Value.Select(i =>
        //    //        (i.MonoisotopicMass + ChargeStateIdentifier.GetDiffToMonoisotopic(i.MassIndex)).ToMz(i.Charge))),
        //    //    string.Join(",", j.Value.Select(i => i.)));
        //});
        //plottingTuples.OrderByDescending(i => i.Item1)
        //    .ForEach(i =>
        //    {
        //        Console.WriteLine(string.Join("\t", i.ToString().Split(separators))); 
        //    });

        //// get the charge states of the isotopic envelopes
        ////var dictChargeStateEnvelopes = GetEnvelopes(results)
        ////    // possible criteria for removing: 
        ////    // non-consecutive charge states; sum of peaks explain less than a certain threshold; 
        ////    // peaks have less
        ////    .OrderByDescending(i => i.Value.Sum(j => j.Score))
        ////    .Take(30)
        ////    .ToDictionary(i => i.Key, i => i.Value); 
        ////// to plot: 
        //// x values: mz of most intense peak in isotopic envelopes 
        //// y values: most intense peak in isotopic envelope 
        //// charge state 
        //// monoisotopic mass 

        //// how do you figure out where to divide the line bewteen real and harmonic? 

        //// plot these three things on top of the original spectrum 
        ////var plottingPoints = GetMzChargeStateIntensityMonoisotopicMass(dictChargeStateEnvelopes);
        ////foreach (var point in plottingPoints)
        ////{
        ////    Console.WriteLine("{0},{1},{2},{3}", point.xVal, point.yVal, point.chargeState, point.monoisotopicMass);
        ////}

    }

    [Test]
    public void TestGroundTruthDataSet()
    {
        string path = Path.Combine(@"C:\Users\Austin\Documents\Projects\MsDataSimulatorOutput\train.mzML");
        MsDataFile file = MsDataFileReader.GetDataFile(path); 
        file.InitiateDynamicConnection();
        MsDataScan scan = file.GetOneBasedScanFromDynamicConnection(116);
        ChargeStateDeconvolutionParams deconParams = new ChargeStateDeconvolutionParams(5, 120, 1, maxThreads: 12);
        ChargeStateIdentifier decon = new ChargeStateIdentifier(deconParams);
        var results = decon.Deconvolute(scan.MassSpectrum, new MzRange(900d, 1000d))
            .ToList();
        var chargeStateEnvelopes = results
            .ObserveAdjacentChargeStates()
            .OrderByDescending(i => i.Key);
        
        // you could score by euclidean distance and then cluster the peaks together i think. 
    }

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
}

//public static class ListIsotopicEnvelopeExtensions
//{
//    public static IEnumerable<IsotopicEnvelope> CleanUpEnvelopes(this IEnumerable<IsotopicEnvelope> envelopes)
//    {
//        foreach (var envelope in envelopes)
//        {
//            if (envelope.Score > 0)
//            {
//                yield return envelope;
//            }
//        }
//    }
//}