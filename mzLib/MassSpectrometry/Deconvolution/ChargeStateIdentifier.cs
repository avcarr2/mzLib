#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using MathNet.Numerics;
using MzLibUtil;
using Constants = Chemistry.Constants;
using Chemistry;
using MathNet.Numerics.Statistics;
using System.Runtime.CompilerServices;

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
        return DeconvolutePrivateFast(spectrumToDeconvolute, range, 0.75);
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
    
    /// <summary>
    /// Given a spectrum and an mz range, uses the ChargeStateLadderMatch object to get the theoretical isotopic envelope, then
    /// grabs the peaks that are within the range peaks in the theoretical isotopic envelope greater than a relative intensity of 0.05. 
    /// </summary>
    /// <param name="match"></param>
    /// <param name="scan"></param>
    /// <param name="range"></param>
    /// <returns></returns>
    internal IEnumerable<IsotopicEnvelope> FindIsotopicEnvelopes(ChargeStateLadderMatch match, MzSpectrum scan, MzRange range)
    {
        this.spectrum = scan;

        for (int i = 0; i < match.ChargesOfMatchingPeaks.Count; i++)
        {
            int charge = (int)Math.Round(match.ChargesOfMatchingPeaks[i]);
            if (range.Contains(match.MatchingMzPeaks[i]))
            {
                double maxMzToTake = allMasses[match.MassIndex].Max().ToMz(charge);
                double minMzToTake = allMasses[match.MassIndex].Min().ToMz(charge);

                MzRange newRange = new(minMzToTake, maxMzToTake);

                var envelope = FillIsotopicEnvelopeByBounds(match, scan, newRange, charge);
                if (envelope != null) yield return envelope; 
            }
        }
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
    internal IEnumerable<ChargeStateLadder> CreateChargeStateLadders(int indexOfMaxIntensityPeak,
        double[] mzValueArray, int minCharge, int maxCharge, double minMzValAllowed, double maxMzValAllowed)
    {
        double mzVal = mzValueArray[indexOfMaxIntensityPeak];
        for (int i = minCharge; i <= maxCharge; i++)
        {
            double tempMass = mzVal.ToMass(i);
            List<(int charge, double mz)> tempLadder = new();
            for (int j = maxCharge; j > 0; j--)
            {
                tempLadder.Add((j, tempMass.ToMz(j)));
            }
            // clean up the ladder 
            tempLadder = tempLadder.Where(k => k.mz <= maxMzValAllowed && k.mz >= minMzValAllowed).ToList();
            yield return new ChargeStateLadder(tempMass, tempLadder.Select(k => k.mz).ToArray());
        }
    }
    /// <summary>
    /// Produces a list of indices that match between a theoretical charge state ladder and the experimental data. 
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="chargeStateLadders"></param>
    /// <param name="ppmMatchTolerance"></param>
    /// <param name="ladderToIndicesMaps"></param>
    private void MatchChargeStateLaddersFast(MzSpectrum scan, List<ChargeStateLadder> chargeStateLadders,
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
    /// <summary>
    /// The function that actually does the work in deconvolution, optimized for speed. 
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="deconvolutionRange"></param>
    /// <param name="spectralSimMatchThresh"></param>
    /// <returns></returns>
    private IEnumerable<IsotopicEnvelope> DeconvolutePrivateFast(MzSpectrum scan, MzRange deconvolutionRange, double spectralSimMatchThresh)
    {
        List<(int index, double mz, double intensity)> mzIntensityPairs =
            (from i in Enumerable.Range(0, scan.XArray.Length)
             select (i, scan.XArray[i], scan.YArray[i]))
            .OrderByDescending(i => i.Item3)
            .ToList();
        List<(int Index, double mz, double intensity)> filteredPairs =
            (from i in Enumerable.Range(0, scan.XArray.Length)
             where scan.XArray[i] >= deconvolutionRange.Minimum && scan.XArray[i] <= deconvolutionRange.Maximum
             select (i, scan.XArray[i], scan.YArray[i]))
            .OrderByDescending(i => i.Item3)
            .ToList();

        // go to next if it is found in seen 

        ConcurrentBag<IsotopicEnvelope> envelopes = new(); 

        Parallel.ForEach(Enumerable.Range(0, filteredPairs.Count), (indexer) =>
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

            foreach (var match in matchList)
            {
                var isotopicEnvelopes = FindIsotopicEnvelopes(match!, scan, deconvolutionRange);
                foreach (var envelope in isotopicEnvelopes)
                {

                    RescoreIsotopicEnvelope(envelope);
                    envelopes.Add(envelope);
                    if (envelope.Score >= spectralSimMatchThresh)
                    {
                        envelopes.Add(envelope);
                    }
                }
            }
        });
        List<IsotopicEnvelope> results = envelopes.OrderByDescending(i => i.Score).ToList();
        
        return results;
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
            .FilterByY(0.01, 1.0);
        var spectrum0 = new MzSpectrum(theoreticalSpectrum.Select(i => i.Mz).ToArray(),
            theoreticalSpectrum.Select(i => i.Intensity).ToArray(), true); 

        MzSpectrum experimentalSpectrum = new MzSpectrum(envelope.Peaks.Select(i => i.mz.ToMass(envelope.Charge)).ToArray(),
            envelope.Peaks.Select(i => i.intensity).ToArray(), true);

        SpectralSimilarity similarity = new SpectralSimilarity(experimentalSpectrum, spectrum0,
            SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak, 1,
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
            // ensures that the lowest mass returned is the monoisotopic mass. 
            //if (envelope.MonoisotopicMass < allMasses[envelope.MassIndex][0] - diffToMonoisotopic[envelope.MassIndex])
            //{
            //    continue;
            //}
            envelope.ScoreByIntensityExplained(scan);
            if(envelope != null) yield return envelope;
        }
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
    
}




