using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using OpenMcdf.Extensions;
using Readers; 
namespace Test;
using NUnit.Framework; 
using MassSpectrometry;

public class TestChargeStateIdentifier
{
    //[Test]
    //public void TestChargeStateIdentification()
    //{
    //    string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
    //        "ExampleChargeStateIdentifierTest.raw");

    //    var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList().First().MassSpectrum;
    //    //scan = ChargeStateIdentifier.DownsampleData(scan, 5); 
    //    var localMaxEnum = ChargeStateIdentifier.ConvertToIndexedTuple(scan.YArray);
    //    var orderedEnum = localMaxEnum.OrderByDescending(i => i.Item2).ToList(); 
    //    var chargeStateLadders = ChargeStateIdentifier.CreateChargeStateLadders(orderedEnum, scan.XArray, 5, 100, 
    //        500, 2000).ToList();
    //    var matchedLadders = ChargeStateIdentifier.MatchChargeStateLadders(scan, chargeStateLadders, orderedEnum, 5, out var ladderToIndicesMap);
    //    var output = ChargeStateIdentifier.TransformToChargeStateLadderMatch(ladderToIndicesMap, matchedLadders, chargeStateLadders.ToList()).ToList();

    //    int indexer = 0; 
    //    // calculate fraction of masses occurring above 0.05 relative intensity
    //    // you're going to have junk values that are additive to the harmonics, but not the true matching value.
    //    // so I'm multiplying the percent of the theoretical peaks identified in the charge state envelope times the intensity explained to
    //    // reduce the intensity explained by the harmonics while preserving that of the original peak. 

    //    var bestScoringChargeStateLadderMatch = ScoreChargeStateLadderMatches(output, scan);
    //    // convert to isotopic envelopes. 
    //        // match to averagine? 

    //    ChargeStateIdentifier csi = new ChargeStateIdentifier(new ClassicDeconvolutionParameters(5, 100, 4, 3));
    //    var result = csi.FindIsotopicEnvelopes(bestScoringChargeStateLadderMatch, scan).ToList();

    //    HashSet<double> usedMzValues = new();
    //    foreach (var envelope in result)
    //    {
    //        if (envelope.Peaks.Count < 5) continue; 
    //        foreach (var mz in envelope.Peaks.Select(i => i.mz))
    //        {
    //            usedMzValues.Add(mz); 
    //        }
    //    }

    //    var intersectedIndices =
    //        (from x in Enumerable.Range(0, scan.XArray.Length)
    //        from y in usedMzValues
    //        where scan.XArray[x] == y
    //        select x)
    //        .OrderBy(i => i)
    //        .ToArray();
    //    List<double> newX = new();
    //    List<double> newY = new();

    //    int intersectionIndex = 0; 
    //    for (int i = 0; i < scan.XArray.Length; i++)
    //    {
    //        if (intersectionIndex < intersectedIndices.Length && intersectedIndices[intersectionIndex] == i)
    //        {
    //            intersectionIndex++;
    //            continue; 
    //        }
    //        else
    //        {
    //            newX.Add(scan.XArray[i]);
    //            newY.Add(scan.YArray[i]);
    //        }
    //    }

    //    var newSpectrum = new MzSpectrum(newX.ToArray(), newY.ToArray(), true); 


    //}
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
        results
            .CleanUpEnvelopes()
            .OrderByDescending(i => i.TotalIntensity)
            //.DistinctBy(i => i.MonoisotopicMass)
            .ForEach(i => 
            { 
                if (i != null)
                {
                    Console.WriteLine("{0}\t{1}\t{2}", i.MonoisotopicMass, i.Charge, i.Score);
                }
            });

    }


    [Test]
    public void TestCreateChargeStateLadder()
    {
        double mass = 1000;
        int charge = 1;
        double adductMass = 1.0078; 

        double exptlMass = (mass + charge * adductMass) / charge;
        double exptlMass2 = mass / charge + adductMass;
        Console.WriteLine(exptlMass + "; " + exptlMass2);
        Console.WriteLine("{0}", (exptlMass - adductMass) * charge);
    }
}

public static class ListIsotopicEnvelopeExtensions
{
    public static IEnumerable<IsotopicEnvelope> CleanUpEnvelopes(this IEnumerable<IsotopicEnvelope> envelopes)
    {
        foreach (var envelope in envelopes)
        {
            if (envelope.Score > 0)
            {
                yield return envelope;
            }
        }
    }
}