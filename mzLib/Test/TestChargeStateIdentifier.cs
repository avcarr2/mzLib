using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics;
using MathNet.Numerics.Statistics;
using Readers; 
namespace Test;
using NUnit.Framework; 
using MassSpectrometry;

public class TestChargeStateIdentifier
{
    [Test]
    public void TestChargeStateIdentification()
    {
        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
            "ExampleChargeStateIdentifierTest.raw");

        var scan = MsDataFileReader.GetDataFile(path).LoadAllStaticData().GetAllScansList().First().MassSpectrum;
        //scan = ChargeStateIdentifier.DownsampleData(scan, 5); 
        var localMaxEnum = ChargeStateIdentifier.ConvertToIndexedTuple(scan.YArray);
        var orderedEnum = localMaxEnum.OrderByDescending(i => i.Item2).ToList(); 
        var chargeStateLadders = ChargeStateIdentifier.CreateChargeStateLadders(orderedEnum, scan.XArray, 5, 100, 
            500, 2000).ToList();
        var matchedLadders = ChargeStateIdentifier.MatchChargeStateLadders(scan, chargeStateLadders, orderedEnum, 5, out var ladderToIndicesMap);
        var output = ChargeStateIdentifier.TransformToChargeStateLadderMatch(ladderToIndicesMap, matchedLadders, chargeStateLadders.ToList()).ToList();

        int indexer = 0; 
        // calculate fraction of masses occurring above 0.05 relative intensity
        // you're going to have junk values that are additive to the harmonics, but not the true matching value.
        // so I'm multiplying the percent of the theoretical peaks identified in the charge state envelope times the intensity explained to
        // reduce the intensity explained by the harmonics while preserving that of the original peak. 

        var bestScoringChargeStateLadderMatch = ScoreChargeStateLadderMatches(output, scan);
        // convert to isotopic envelopes. 
            // match to averagine? 
        


    }

    public ChargeStateLadderMatch ScoreChargeStateLadderMatches(List<ChargeStateLadderMatch> ladderMatches, MzSpectrum scan)
    {
        return ladderMatches.OrderByDescending(i =>
            ScoreByIntensityExplained(i, scan) * CompareTheoreticalNumberChargeStatesVsActual(i)).First();
    }

    public void CleanUpChargeStateLadderMatch(ChargeStateLadderMatch match, MzSpectrum spectra)
    {
        // examine the end and remove if non-consecutive. 
        double[] diffs = new double[match.ChargesOfMatchingPeaks.Count - 1];
        // if diff is greater than 1 => you need to check to see if the other, intermediate charge states exist
        // criteria for if a charge state 
        for (int i = 0; i < diffs.Length; i++)
        {
            double diff = match.ChargesOfMatchingPeaks[i] - match.ChargesOfMatchingPeaks[i + 1]; 

            double tolerance = Math.Abs(1 - Math.Round(diff));
            if (tolerance < 0.1)
            {
                continue;
            }else if (tolerance >= 2)
            {
                match.ChargesOfMatchingPeaks.RemoveAt(i + 1);
                match.IntensitiesOfMatchingPeaks.RemoveAt(i + 1);
                match.MatchingMzPeaks.RemoveAt(i + 1);
            }

            double chargeState1 = match.ChargesOfMatchingPeaks[i];
            double chargeState2 = match.ChargesOfMatchingPeaks[i + 1];
            while (chargeState1 > chargeState2)
            {
                // look for the intermediate peaks in the original data. 
                //Array.BinarySearch()
                
                //spectra.XArray

                //chargeState1--; 
            }

        }


        // examine the middle and see if there are any missing charge states. Check the data at an increased tolerance to see if they exist. 
        
    }

    public void QualityControlOutputEnvelopes(List<ChargeStateLadderMatch> ladderMatches)
    {

    }

    public double CompareTheoreticalNumberChargeStatesVsActual(ChargeStateLadderMatch ladderMatch)
    {
        return ladderMatch.ChargesOfMatchingPeaks.Count / (double)ladderMatch.TheoreticalLadder.MzVals.Length; 
    }

    public double ScoreSequentialChargeStates(ChargeStateLadderMatch ladderMatch)
    {
        double[] output = new double[ladderMatch.ChargesOfMatchingPeaks.Count - 1]; 
        for (int i = 0; i < ladderMatch.ChargesOfMatchingPeaks.Count - 1; i++)
        {
            output[i] = ladderMatch.ChargesOfMatchingPeaks[i] - ladderMatch.ChargesOfMatchingPeaks[i + 1]; 
        }

        return 1 / output.Sum() / (double)output.Length; 
    }


    public double ScoreByIntensityExplained(ChargeStateLadderMatch match, MzSpectrum spectrum, double threshold = 0.05)
    {
        return match.IntensitiesOfMatchingPeaks.Select(i => i / spectrum.YArray.Max()).Where(i => i >= threshold).Sum();
    }

    public void CalculateMasses(ChargeStateLadderMatch match)
    {
        double[] masses = new double[match.ChargesOfMatchingPeaks.Count];
        for (int i = 0; i < match.ChargesOfMatchingPeaks.Count; i++)
        {
            //match.ChargesOfMatchingPeaks[i]
        }
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