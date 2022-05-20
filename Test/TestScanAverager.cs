using InstrumentControl;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    public static class TestScanAverager
    {

		// TESTING IDEA: take scans, group them into five, send them through with different parametes
		// and run a classic search on teh resulting composit scans. Print a txt file that has the 
		// RejectionType, WeightingType, (paramertes of weighting if applicable), PSMs, peptides, proteins
		// also save the psmtsv from the searches with identifying names.
		// This should allow a simple way to optimize parameters

		// Used to initialte the main method of the class, will need to turn into an actual test in the future
		[Test]
		public static void TestSpectrumBinning()
		{
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			List<MsDataScan> scans = MS1DatabaseParser.LoadAllScansFromFile(filepath);

			
			var reducedScans = scans.Where(p => p.MsnOrder == 1).TakeLast(5).ToList();
			ScanAverager.AverageScans(reducedScans);



			ScanAverager.ResetValues();
		}

		[Test]
		public static void TestClippingFunctions()
        {
			#region Min Max Clipping

			double[] test = { 10, 9, 8, 7, 6, 5 };
			double[] expected = { 9, 8, 7, 6 };
			double[] minMaxClipped = ScanAverager.MinMaxClipping(test);
			Assert.AreEqual(expected, minMaxClipped);

			#endregion

			#region Percentile Clipping

			test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			expected = new double[] { 90, 80, 70, 60, 50, 40, 30, 20, 10 };
			double[] percentileClipped = ScanAverager.PercentileClipping(test, 0.9);
			Assert.AreEqual(expected, percentileClipped);

			test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
			expected = new double[] { 60, 50, 40, 30, 20, 10 };
			percentileClipped = ScanAverager.PercentileClipping(test, 0.9);
			Assert.AreEqual(expected, percentileClipped);

			#endregion

			#region Sigma Clipping

			test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
			double[] sigmaClipped = ScanAverager.SigmaClipping(test, 1.5, 1.5);
			expected = new double[] { 50, 40, 30, 20, 10 };
			Assert.AreEqual(expected, sigmaClipped);

			test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			sigmaClipped = ScanAverager.SigmaClipping(test, 1, 1);
			expected = new double[] { 50 };
			Assert.AreEqual(expected, sigmaClipped);

			sigmaClipped = ScanAverager.SigmaClipping(test, 1.3, 1.3);
			expected = new double[] { 60, 50, 40 };
			Assert.AreEqual(expected, sigmaClipped);

			sigmaClipped = ScanAverager.SigmaClipping(test, 1.5, 1.5);
			expected = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			Assert.AreEqual(expected, sigmaClipped);

			#endregion

			#region Windsorized Sigma Clipping

			test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
			double[] windsorizedSigmaClipped = ScanAverager.WinsorizedSigmaClipping(test, 1.5, 1.5);
			expected = new double[] { 60, 50, 40, 30, 20, 10, 0 };
			Assert.AreEqual(expected, windsorizedSigmaClipped);

			test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			windsorizedSigmaClipped = ScanAverager.WinsorizedSigmaClipping(test, 1, 1);
			expected = new double[] { 60, 50, 40 };
			Assert.AreEqual(expected, windsorizedSigmaClipped);

			windsorizedSigmaClipped = ScanAverager.WinsorizedSigmaClipping(test, 1.3, 1.3);
			expected = new double[] { 60, 50, 40 };
			Assert.AreEqual(expected, windsorizedSigmaClipped);

			windsorizedSigmaClipped = ScanAverager.WinsorizedSigmaClipping(test, 1.5, 1.5);
			expected = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			Assert.AreEqual(expected, windsorizedSigmaClipped);

			#endregion

			#region Averaged Sigma Clipping

			test = new double[] { 100, 80, 60, 50, 40, 30, 20, 10, 0 };
			double[] averagedSigmaClipping = ScanAverager.AveragedSigmaClipping(test, 3, 3);
			expected = new double[] { 80, 60, 50, 40, 30, 20, 10, 0 };
			Assert.AreEqual(expected, averagedSigmaClipping);

			test = new double[] { 100, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0 };
			averagedSigmaClipping = ScanAverager.AveragedSigmaClipping(test, 1, 1);
			expected = new double[] { 70, 60, 50, 40, 30, };
			Assert.AreEqual(expected, averagedSigmaClipping);

			averagedSigmaClipping = ScanAverager.AveragedSigmaClipping(test, 1.3, 1.3);
			expected = new double[] { 80, 70, 60, 50, 40, 30, 20 };
			Assert.AreEqual(expected, averagedSigmaClipping);

			averagedSigmaClipping = ScanAverager.AveragedSigmaClipping(test, 1.5, 1.5);
			expected = new double[] { 80, 70, 60, 50, 40, 30, 20 };
			Assert.AreEqual(expected, averagedSigmaClipping);

			#endregion

		}

		[Test]
		public static void TestWeightedAverage()
        {
			double[] values = new double[] { 10, 0 };
			double[] weights = new double[] { 8, 2 };
			double average = ScanAverager.CalculateAverage(values, weights);
			Assert.AreEqual(8, average);

			values = new double[] { 10, 2, 0 };
			weights = new double[] { 9, 1, 0 };
			average = ScanAverager.CalculateAverage(values, weights);
			Assert.AreEqual(9.200, Math.Round(average,4));
        }

		[Test]
		public static void TestCalculatingWeights()
        {
			double[] test = new double[] { 10, 8, 6, 5, 4, 3, 2, 1 };
			double[] weights = new double[test.Length];
			ScanAverager.WeightByNormalDistribution(test, ref weights);
			double weightedAverage = ScanAverager.CalculateAverage(test, weights);
			Assert.AreEqual(4.5460 , Math.Round(weightedAverage,4));

			weights = new double[test.Length];
			ScanAverager.WeightByCauchyDistribution(test, ref weights);
			weightedAverage = ScanAverager.CalculateAverage(test, weights);
			Assert.AreEqual(4.6411, Math.Round(weightedAverage, 4));

			weights = new double[test.Length];
			ScanAverager.WeightByPoissonDistribution(test, ref weights);
			weightedAverage = ScanAverager.CalculateAverage(test, weights);
			Assert.AreEqual(4.2679, Math.Round(weightedAverage, 4));

			weights = new double[test.Length];
			ScanAverager.WeightByGammaDistribution(test, ref weights);
			weightedAverage = ScanAverager.CalculateAverage(test, weights);
			Assert.AreEqual(4.1638, Math.Round(weightedAverage, 4));
        }

		[Test]
		public static void TestSetValuesAndRejectOutliersSwitch()
        {
			ScanAverager.ResetValues();
			Assert.That(ScanAverager.RejectionType == RejectionType.NoRejection);
			Assert.AreEqual(ScanAverager.Percentile, 0.9);
			Assert.AreEqual(ScanAverager.MinSigmaValue, 1.3);
			Assert.AreEqual(ScanAverager.MaxSigmaValue, 1.3);

			ScanAverager.SetValues(RejectionType.MinMaxClipping, WeightingType.NoWeight, .8, 2, 4);
			Assert.That(ScanAverager.RejectionType == RejectionType.MinMaxClipping);
			Assert.AreEqual(ScanAverager.Percentile, 0.8);
			Assert.AreEqual(ScanAverager.MinSigmaValue, 2);
			Assert.AreEqual(ScanAverager.MaxSigmaValue, 4);

			ScanAverager.ResetValues();
			double[] test = new double[] { 10, 8, 6, 5, 4, 2, 0 };
			double[] output = ScanAverager.RejectOutliers(test);
			double[] expected = new double[] { 10, 8, 6, 5, 4, 2, 0 };
			Assert.AreEqual(expected, output);

			ScanAverager.SetValues(rejectionType: RejectionType.MinMaxClipping);
			Assert.That(ScanAverager.RejectionType == RejectionType.MinMaxClipping);
			output = ScanAverager.RejectOutliers(test);
			expected = new double[] { 8, 6, 5, 4, 2 };
			Assert.AreEqual(expected, output);

			ScanAverager.SetValues(rejectionType: RejectionType.PercentileClipping);
			Assert.That(ScanAverager.RejectionType == RejectionType.PercentileClipping);
			ScanAverager.SetValues(rejectionType: RejectionType.SigmaClipping);
			Assert.That(ScanAverager.RejectionType == RejectionType.SigmaClipping);
			ScanAverager.SetValues(rejectionType: RejectionType.WinsorizedSigmaClipping);
			Assert.That(ScanAverager.RejectionType == RejectionType.WinsorizedSigmaClipping);
			ScanAverager.SetValues(rejectionType: RejectionType.AveragedSigmaClipping);
			Assert.That(ScanAverager.RejectionType == RejectionType.AveragedSigmaClipping);
			ScanAverager.ResetValues();
		}

	}
}
