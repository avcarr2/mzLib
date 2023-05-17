using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using OxyPlot.Wpf;
using Plotly.NET;
using Proteomics.AminoAcidPolymer;
using SimulatedData;
using UsefulProteomicsDatabases;
// remove plotly before final merge 
using Plotly.NET.CSharp;
using Chart = Plotly.NET.CSharp.Chart;

namespace Test
{
	public class TestSimulatedChargeStateEnvelope
	{
		// phospholamban sequence, UniProt entry P26678
		private const string actin = "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEKIWHHTFYNELRVAPEEHPVLLTEAPLNPKANREKMTQIMFETFNTPAMYVAIQAVLSLYASGRTTGIVMDSGDGVTHTVPIYEGYALPHAILRLDLAGRDLTDYLMKILTERGYSFTTTAEREIVRDIKEKLCYVALDFEQEMATAASSSSLEKSYELPDGQVITIGNERFRCPEALFQPSFLGMESCGIHETTFNSIMKCDVDIRKDLYANTVLSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWISKQEYDESGPSIVHRKCF";

		private const string cytochromeC =
			"MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"; 
		private Chemistry.IsotopicDistribution _distribution1;
		private IsotopicDistribution _distribution2; 
		private SimulatedChargeStateEnvelope Cse;
		static double mzLow = 500d;
		static double mzHigh = 2000d;
		static double stepSize = 0.001;
		static int length = (int)((mzHigh - mzLow) / stepSize);
		static int chargeStateLow = 5;
		static int chargeStateHigh = 100;
		static (double, double) chargeStateDistr = (-0.65, 0.15);

		[OneTimeSetUp]
		public void OneTimeSetup()
		{
			UsefulProteomicsDatabases.Loaders.LoadElements(); ; 
			_distribution1 = IsotopicDistribution.GetDistribution(new Peptide(actin).GetChemicalFormula());
			_distribution2 = IsotopicDistribution.GetDistribution(new Peptide(cytochromeC).GetChemicalFormula());
		}

		[Test]
		public void TestIncorrectInstatiation()
		{
			
			// checks that error handling of initialization accurately throws 
			// an mzLib exception. 
			Assert.Throws<MzLibException>(() =>
			{
				new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length,
					chargeStateHigh, chargeStateLow, 
					_distribution1, chargeStateDistr);
			}); 
			Assert.Throws<MzLibException>(() =>
			{
				new SimulatedChargeStateEnvelope(mzHigh, mzLow, stepSize, length,
					chargeStateLow, chargeStateHigh, 
					_distribution1, chargeStateDistr);
			});
		}

		[Test]
		public void TestSimulatedChargeStateInstantiation()
		{
			Cse = new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length,
				chargeStateLow, chargeStateHigh,
				_distribution1, chargeStateDistr);
			//Chart.Line<double, double, string>(Cse.Xarray, Cse.Yarray).Show();
			// expected values
			List<(double, double)> expectedList = new List<(double, double)>
			{
				// max peak in the spectra
				(1222.487, 0.78739), 
				// second largest peak
				(1222.287, 0.5703315)
			};
			double expectedMzDifference = expectedList[0].Item1 - expectedList[1].Item1;
			int expectedChargeState = (int)Math.Round(1d / expectedMzDifference);

			// values from instantiation
			double[] top2IntensityVals = Cse.Yarray.OrderByDescending(i => i).Take(2).ToArray();
			int[] top2IntIndex = top2IntensityVals.Select(i => (int)Cse.Yarray.IndexOf(i)).ToArray();
			double[] top2MzVals = new double[] { Cse.Xarray[top2IntIndex[0]], Cse.Xarray[top2IntIndex[1]] };
			int actualChargeState = (int)Math.Round(1d / (top2MzVals[0] - top2MzVals[1]));

			Assert.AreEqual(expectedChargeState, actualChargeState);
			Assert.That(new[]{ expectedList[0].Item2, expectedList[1].Item2 }, Is.EqualTo(top2IntensityVals).Within(0.001));
		}
		[Test]
		[Repeat(5)]
		public void TestAddHighFrequencyNoise()
		{
			Cse = new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length,
				chargeStateLow, chargeStateHigh,
				_distribution1, chargeStateDistr);
			// hf noise distribution
			Normal hfNoiseDistr = new(0.1, 0.01); 
			
			Cse.AddHighFrequencyNoise(hfNoiseDistr);
			// check that the first 100 values have a mean of 0.1 and a std 0f 0.01 
			(double, double) meanStddev = Cse.Yarray[..100].MeanStandardDeviation(); 
			Assert.That(meanStddev, Is.EqualTo((0.1, 0.01)).Within(0.01));
		}

		[Test]
		[Repeat(5)]
		public void TestAddLowFrequencyNoise()
		{
			LowFrequencyNoiseParameters lfNoiseParams = new(5, 10, 501d, 1999d,
				0.5, 1, 0.01, 0.2);
			Cse = new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length,
				chargeStateLow, chargeStateHigh,
				_distribution1, chargeStateDistr);

			double[] originalYarrayCopy = new double[Cse.Yarray.Length];
			Buffer.BlockCopy(Cse.Yarray, 0, 
				originalYarrayCopy, 0, sizeof(double) * Cse.Yarray.Length);

			List<int> seedsList = new List<int>()
			{
				1551, 1552, 1553, 1554, 1555
			};
			Cse.AddLowFrequencyNoise(lfNoiseParams, seedsList);
			// subtract the orignal y array from the yarray with lf noise peaks added to it. 
			// check that the maximum noise peak is at the correct location. 
			double[] lfNoiseArray = Cse.Yarray
				.Select((i, index) => i - originalYarrayCopy[index])
				.ToArray();
			double maxLfNoisePeak = lfNoiseArray.Max();
			Assert.That(maxLfNoisePeak, Is.EqualTo(0.1015).Within(0.005));
		}

		[Test]
		public void TestMultipleOverlappingProteoforms()
		{
			var cseCytochromeC = new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length, 5, 25,
				_distribution2, chargeStateDistr);
			var cseActin = new SimulatedChargeStateEnvelope(mzLow, mzHigh, stepSize, length, 25, 60, _distribution1,
				chargeStateDistr);

			//var cytoPlot = Chart.Line<double, double, string>(cseCytochromeC.Xarray, cseCytochromeC.Yarray,
			//	Name: "UniProt Accession: P60709"); 
			//var actinPlot = Chart.Line<double, double, string>(cseActin.Xarray, cseActin.Yarray, Name:"UniProt Accession: P60709");

			//var combined = new List<Plotly.NET.GenericChart.GenericChart>()
			//{
			//	cytoPlot,
			//	actinPlot
			//};  

			//Plotly.NET.GenericChartExtensions.Show(GenericChart.combine(combined));

			using var cytoWriter = new StreamWriter("cytoSimData.txt");
			using var actinWriter = new StreamWriter("actinSimulated.txt");

			for (int i = 0; i < cseCytochromeC.Xarray.Length; i++)
			{
				cytoWriter.WriteLine(cseCytochromeC.Xarray[i] + "\t" + cseCytochromeC.Yarray[i]);
			}
			cytoWriter.Flush();
			for (int i = 0; i < cseActin.Xarray.Length; i++)
			{
				actinWriter.WriteLine(cseActin.Xarray[i] + "\t" + cseActin.Yarray[i]);
			}
			actinWriter.Flush();
		}

		[Test]
		public void TestNormalizeByTic()
		{

		}

		[Test]
		public void TestToMzSpectrum()
		{

		}

		[Test]
		public void TestToMsDataScan()
		{

		}

		[Test]
		public void TestBlur()
		{

		}

		[Test]
		public void NormalizeToMaxIntensity()
		{

		}
	}
}
