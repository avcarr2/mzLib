using MathNet.Numerics.Statistics;
using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class ChargeStateLadderMatch
    {
        public ChargeStateLadder TheoreticalLadder { get; set; }
        public List<double> MatchingMzPeaks { get; set; }
        public List<double> IntensitiesOfMatchingPeaks { get; set; }
        public List<double> ChargesOfMatchingPeaks { get; set; }
        public double EnvelopeScore { get; set; }
        public double Score { get; set; }
        public int MassIndex { get; set; }
        public double MonoisotopicMass => MonoGuesses.Average();
        public double StdDev => MonoGuesses.StandardDeviation(); 
        public double SequentialChargeStateScore { get; set; }
        public List<double> MonoGuesses { get; set; } = new();
        public List<double> MonoErrors { get; set; } = new();
        public double PercentageMzValsMatched { get; set; }

        internal void CalculateMatchError()
        {
            if (!MatchingMzPeaks.Any()) return;

            var orderedMzVal = TheoreticalLadder.MzVals.OrderBy(x => x).ToArray();
            foreach (var peak in MatchingMzPeaks)
            {
                int index = ChargeStateIdentifier.GetBucket(orderedMzVal, peak);
                double ppmError = Math.Abs((peak - orderedMzVal[index]) / peak * 1E6);
                MonoErrors.Add(ppmError);
            }
        }
        

        internal void RemoveHighErrorValues(double errorThreshold)
        {
            for (int i = 0; i < MonoErrors.Count; i++)
            {
                if (MonoErrors[i] > errorThreshold)
                {
                    // remove at monoErrors, Intensities, Charges
                    // note that monoGuesses gets populated later
                    MonoErrors.RemoveAt(i);
                    MatchingMzPeaks.RemoveAt(i);
                    IntensitiesOfMatchingPeaks.RemoveAt(i);
                    ChargesOfMatchingPeaks.RemoveAt(i);
                    i--; 
                }
            }
        }

        /// <summary>
        /// Calculates the total fraction of intensity explained by the mz values that match to the theoretical envelope. 
        /// </summary>
        /// <param name="spectrum"></param>
        /// <param name="threshold"></param>
        /// <returns></returns>
        internal void ScoreByIntensityExplained(MzSpectrum spectrum, double threshold = 0.00001)
        {
            double spectrumThresholdVal = spectrum.YArray.Max() * threshold;
            Score = IntensitiesOfMatchingPeaks.Where(i => i >= spectrumThresholdVal).Sum();
        }
        /// <summary>
        /// Calculates the average spacing between distinct charge states. If the charge states are sequential,
        /// then the charges will be close to -1. If the charge states correspond to a high harmonic, then the values of the charges states are going to be
        /// less than -1.
        /// </summary>
        internal void CalculateChargeStateScore()
        {

            var chargesList = ChargesOfMatchingPeaks
                .Zip(ChargesOfMatchingPeaks.Skip(1), (x, y) => y - x)
                // Where clause removes any diffs that are less than 0.1, which would indicate that the 
                // peak matching tolerance was broad enough to grab more than one peak. However, because the 
                // ppm tolerance for peak matching is so low, the difference between consecutive peaks should remain 
                // approximately the same. 
                .Where(i => Math.Abs(i) > 0.1).ToList();
            if (chargesList.Any())
            {
                var numberUniqueChargesWithCorrectChargeState = chargesList
                    .Select(i => Math.Abs(i - 1d))
                    .Count(i => i <= 0.1); 

                SequentialChargeStateScore = (double)numberUniqueChargesWithCorrectChargeState / (double)chargesList.Count;
                return;
            }

            SequentialChargeStateScore = -10000;
        }
        /// <summary>
        /// The envelope score determines the coefficient of best fit between the (mz,intensity) values of the
        /// matched peaks to a second order polynomial. This will eliminate peaks that are primarily made of noise values and
        /// peaks that don't have enough points to calculate the score. 
        /// </summary>
        internal void CalculateEnvelopeScore()
        {
            if (IntensitiesOfMatchingPeaks.Count < 4)
            {
                EnvelopeScore = 0d;
                return;
            }
            double[] coefficients =
                Fit.Polynomial(ChargesOfMatchingPeaks.ToArray(), IntensitiesOfMatchingPeaks.ToArray(), 2);
            double c = coefficients[0];
            double b = coefficients[1];
            double a = coefficients[2];

            // calculate theoretical polynomial to get R^2. 
            double[] theoreticalPolynom = ChargesOfMatchingPeaks
                .Select(x => c + b * x + a * x * x)
                .ToArray();
            double sum1 = 0;
            double sum2 = 0;
            double ymean = IntensitiesOfMatchingPeaks.Mean();

            for (int i = 0; i < ChargesOfMatchingPeaks.Count; i++)
            {
                sum1 += Math.Pow(IntensitiesOfMatchingPeaks[i] - theoreticalPolynom[i], 2);
                sum2 += Math.Pow(IntensitiesOfMatchingPeaks[i] - ymean, 2);
            }

            EnvelopeScore = 1 - sum1 / sum2;
        }
        /// <summary>
        /// Calculates the percentage of matched mz values in an experimental spectrum compared to a set of theoretical mz values. 
        /// </summary>
        /// <param name="ladderMatch"></param>
        /// <returns></returns>
        internal void CompareTheoreticalNumberChargeStatesVsActual()
        {
            int integerUniqueChargeValuesLength = ChargesOfMatchingPeaks
                .Select(i => (int)Math.Round(i))
                .Distinct()
                .Count();

            PercentageMzValsMatched = (double)integerUniqueChargeValuesLength / (double)TheoreticalLadder.MzVals.Length;
        }
    }

}
