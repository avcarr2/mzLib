﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography.Xml;
using System.Windows.Markup;
using System.Windows.Media.Animation;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Readers; 
namespace Test;
using NUnit.Framework; 
using MassSpectrometry;

public class TestChargeStateIdentifier
{


    [Test]
    [TestCase(2,60)]
    public void TestDeconvolution(int minCharge, int maxCharge)
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "chargeStateDeconTest3raw.raw");
        FilteringParams filteringParams = new FilteringParams(minimumAllowedIntensityRatioToBasePeak:0.01);
        var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData(filteringParams).GetAllScansList().First().MassSpectrum;

        ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 5); 
        ChargeStateIdentifier csi = new(deconParams);

        var results = csi.Deconvolute(scan, scan.Range).ToList();
        results.ToList()
            .OrderByDescending(i => i.TotalIntensity)
            //.DistinctBy(i => i.MonoisotopicMass)
            .ForEach(i => 
            { 
                if (i != null)
                {
                    //Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
                }
            });

        // get the charge states of the isotopic envelopes
        var dictChargeStateEnvelopes = GetEnvelopes(results)
            // possible criteria for removing: 
            // non-consecutive charge states; sum of peaks explain less than a certain threshold; 
            // peaks have less
            .OrderByDescending(i => i.Value.Sum(j => j.Score))
            .Take(30)
            .ToDictionary(i => i.Key, i => i.Value); 
        // to plot: 
        // x values: mz of most intense peak in isotopic envelopes 
        // y values: most intense peak in isotopic envelope 
        // charge state 
        // monoisotopic mass 

        // how do you figure out where to divide the line bewteen real and harmonic? 

        // plot these three things on top of the original spectrum 
        var plottingPoints = GetMzChargeStateIntensityMonoisotopicMass(dictChargeStateEnvelopes);
        foreach (var point in plottingPoints)
        {
            Console.WriteLine("{0},{1},{2},{3}", point.xVal, point.yVal, point.chargeState, point.monoisotopicMass);
        }

    }

    [Test]
    public void TestDeconvolutionBigThings()
    {
        
        int minCharge = 5;
        int maxCharge = 100; 
        string path = @"D:\DeconDataSet\SEC4-08AUG16_5uLinj_3SEC_000021 (2)..mzML";
        var reader = MsDataFileReader.GetDataFile(path); 
        reader.InitiateDynamicConnection();
        var scan = reader.GetOneBasedScanFromDynamicConnection(679).MassSpectrum;

        ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 20);
        ChargeStateIdentifier csi = new(deconParams);

        var results = csi.Deconvolute(scan, new MzRange(600,650)).ToList();
        
        results
           // .CleanUpEnvelopes()
            .OrderByDescending(i => i.Score)
            //.DistinctBy(i => i.MonoisotopicMass)
            .ForEach(i =>
            {
                if (i != null)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
                }
            });

        // get the charge states of the isotopic envelopes
        var dictChargeStateEnvelopes = GetEnvelopes(results)
            // possible criteria for removing: 
            // non-consecutive charge states; sum of peaks explain less than a certain threshold; 
            // peaks have less
            .OrderByDescending(i => i.Value.Sum(j => j.Score))
            .Take(30)
            .ToDictionary(i => i.Key, i => i.Value);
        // to plot: 
    }

    [Test]
    public void VeryBigProteinDeconvolutionTest()
    {
        int minCharge = 50;
        int maxCharge = 250;
        string path = @"D:\DeconDataSet\SEC4-08AUG16_5uLinj_3SEC_000021 (2)..mzML";
        var reader = MsDataFileReader.GetDataFile(path);
        reader.InitiateDynamicConnection();
        var scan = reader.GetOneBasedScanFromDynamicConnection(1016).MassSpectrum;

        ChargeStateDeconvolutionParams deconParams = new(minCharge, maxCharge, 10);
        ChargeStateIdentifier csi = new(deconParams);

        var results = csi.Deconvolute(scan, new MzRange(985, 1000)).ToList();

        results
            // .CleanUpEnvelopes()
            .OrderByDescending(i => i.TotalIntensity)
            //.DistinctBy(i => i.MonoisotopicMass)
            .ForEach(i =>
            {
                if (i != null)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
                }
            });
        // suppose that you have isolated a 

    }

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
                    (ChargeStateIdentifier.GetDiffToMonoisotopic(envelope.MassIndex) + mass) / (double)envelope.Charge + Chemistry.Constants.ProtonMass;
                // retrieve the intensity corresponding to the mostabundant isotopic mass
                double[] arrayOfMzVals = envelope.Peaks
                    .OrderBy(i => i.mz)
                    .Select(i => i.mz)
                    .ToArray(); 

                int indexOfMostAbundant = Array.BinarySearch(arrayOfMzVals, mostAbundantMz); 

                indexOfMostAbundant = indexOfMostAbundant < 0 ? ~indexOfMostAbundant : indexOfMostAbundant;
                indexOfMostAbundant = indexOfMostAbundant >= envelope.Peaks.Count ? envelope.Peaks.Count - 1 : indexOfMostAbundant;
                // the OrderBy shouldn't be repeated, but I'm tired and about to go home. 
                var maxPairs = envelope.Peaks.OrderBy(i=>i.mz).ToArray()[indexOfMostAbundant]; 

                int charge = envelope.Charge;
                outputList.Add((maxPairs.mz, maxPairs.intensity, charge, mass));
            }
        }
        return outputList;
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