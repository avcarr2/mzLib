using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using FtmsSimulator;
using MassSpectrometry;
using Proteomics;
using Chemistry;
using OxyPlot;
using System.IO; 

namespace Test
{
	class TestFTMSSimulator
	{
        public static void WriteOxyPlotToPDF(string fileName, PlotModel plotmodel)
        {
            using (var stream = File.Create(fileName))
            {
                var pdfExporter = new PdfExporter { Width = 600, Height = 600 };
                pdfExporter.Export(plotmodel, stream);
            }
        }
        [Test]
        public void TestSimulatedProteinInitiation()
        {
            // SERCA2a sequence: 
            /*string sequence = "MENAHTKTVEEVLGHFGVNESTGLSLEQVKKLKERWGSNELPAEEGKTLLELVIEQFEDLLVRILLLAACISFVLAWFEEGEET" +
                "ITAFVEPFVILLILVANAIVGVWQERNAENAIEALKEYEPEMGKVYRQDRKSVQRIKAKDIVPGDIVEIAVGDKVPADIRLTSIKSTTLRVDQSILTGE" +
                "SVSVIKHTDPVPDPRAVNQDKKNMLFSGTNIAAGKAMGVVVATGVNTEIGKIRDEMVATEQERTPLQQKLDEFGEQLSKVISLICIAVWIINIGHFNDP" +
                "VHGGSWIRGAIYYFKIAVALAVAAIPEGLPAVITTCLALGTRRMAKKNAIVRSLPSVETLGCTSVICSDKTGTLTTNQMSVCRMFILDRVEGDTCSLNEF" +
                "TITGSTYAPIGEVHKDDKPVNCHQYDGLVELATICALCNDSALDYNEAKGVYEKVGEATETALTCLVEKMNVFDTELKGLSKIERANACNSVIKQLMKKEF" +
                "TLEFSRDRKSMSVYCTPNKPSRTSMSKMFVKGAPEGVIDRCTHIRVGSTKVPMTSGVKQKIMSVIREWGSGSDTLRCLALATHDNPLRREEMHLEDSANFIKY" +
                "ETNLTFVGCVGMLDPPRIEVASSVKLCRQAGIRVIMITGDNKGTAVAICRRIGIFGQDEDVTSKAFTGREFDELNPSAQRDACLNARCFARVEPSHKSKIVEFLQSFDE" +
                "ITAMTGDGVNDAPALKKAEIGIAMGSGTAVAKTASEMVLADDNFSTIVAAVEEGRAIYNNMKQFIRYLISSNVGEVVCIFLTAALGFPEALIPVQLLWVNLVTDGLPATALGFNP" +
                "PDLDIMNKPPRNPKEPLISGWLFFRYLAIGCYVGAATVGAAAWWFIAADGGPRVSFYQLSHFLQCKEDNPDFEGVDCAIFESPYPMTMALSVLVTIEMCNALNSLSENQSLLRMPPWENIWLV" +
                "GSICLSMSLHFLILYVEPLPLIFQITPLNVTQWLMVLKISLPVILMDETLKFVARNYLEPGKECVQPATKSCSFSACTDGISWPFVLLIMPLVIWVYSTDTNFSDMFWS";*/
            
            // Generate the expected isotopic distribution for a fragment of SERCA2a sequence: 
            string sequence = "MENAHTKTVEEVLGHFGVNESTGLSLEQVKKLKERWGSNELPAEEGKTLLELVIEQFEDLLVRILLLAACISFVLAWFEEGEET"; 
            SimulatedProtein simProt = new SimulatedProtein(sequence);
            ChemicalFormula simProtFormula = simProt.GetChemicalFormula();
            Console.WriteLine(simProtFormula.Formula);
            IsotopicDistribution isoEnvelope = IsotopicDistribution.GetDistribution(simProtFormula);
            var plotModel = SimulatedProteinPlot.PlotIsotopicDistribution(isoEnvelope);
            WriteOxyPlotToPDF("isotopic_envelope.pdf", plotModel);

            // get the masses and intensities for each mass in an array. 
                // should probably put this into their own fields within the SimulatedProtein object. Would be easier to access. 
            double[] masses = isoEnvelope.Masses.ToArray();
            double[] intensities = isoEnvelope.Intensities.ToArray();

            // Lines 56-60 need to be in their own method. 
            // for each number of ions, calculate the isotopic envelope for each charge state
            int[] chargeStates = Enumerable.Range(35,11).ToArray();
            // this was a calculated value
            double[] ionCount = { 9766, 43945, 117188, 205078, 246094, 205078, 117188, 43945, 9466, 977, 89 };
            // need to ion count * intensities to get the number of ions per isotope per charge state
                //Wrote method for the below code, but this works and I haven't tested the new method yet.
            List<(int, double, double)> ionCountsTuple = new(); 
            for(int i = 0; i < ionCount.Length; i++)
			{
                for(int j = 0; j < masses.Length; j++)
				{
                    double tempIonCount = intensities[j] * ionCount[i];
                    ionCountsTuple.Add((chargeStates[i], masses[j], tempIonCount)); 
				}
			}

            //Constants: 
            // Calculate the ion decay coefficient for the entire packet of ions: 
            double minT = 0.0;
            double maxT = 2.0;
            double deltaT = 0.000001;
            double percentDecay = 0.6;
            double decayConstant = Transient.CalculateIonDecayConstant(maxT, percentDecay, (int)ionCount.Sum()); 
            // Build list of transients from all elements in the ionCountsTuple list
            List<Transient> transientList = new List<Transient>();
            InstrumentParameters instrPara = new(0.05, 1.0, 300000, true, 0.80);
            Parallel.ForEach(ionCountsTuple, i =>
            {
                var tempTransient = new Transient(i, instrPara, minT, maxT, deltaT, true, percentDecay, decayConstant);
                transientList.Add(tempTransient); 
            });
            double[] summedTransient = Transient.SumMultipleTransients(transientList); 
            Console.WriteLine(string.Join("\n", summedTransient));
            WriteOxyPlotToPDF("transientplot.pdf", TransientPlot.PlotTransient(summedTransient)); 

            // Using too many casts in this code.. I need to update some method signatures to solve that issue... 
        }
        [Test]

        public void TestRange()
		{
            Console.WriteLine(string.Join(",", Enumerable.Range(30, 11).ToArray()));
		}

        [Test]
        public void TestDistribution()
		{
            
            int[] seq = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
            double[] expResults = SimulatedProtein.CalculateBinomialDistribution(seq, 10, 0.5);
            Console.WriteLine(string.Join(", ", expResults)); 
		}
        [Test]
        public void TestBinomialDistribution()
		{
            int[] seq = { 1, 2, 3, 4, 5, 6, 7 };
            int fact1 = SimulatedProtein.Factorial(9);
            int fact2 = SimulatedProtein.Factorial(8);
            int result = fact1 - fact2;
            int result2 = SimulatedProtein.Factorial(10 - 7);
            Console.Write(SimulatedProtein.Factorial(100).ToString()); 
            Console.WriteLine("{0}, {1}", result, result2); 

		}
        [Test]
        public void TestTransientSimulation()
		{
            var testTransient = new Transient(200, 2, 2, 1000000, 10, 1.0, 0, 1, 0.001, false, false, 1.0);
            Console.WriteLine(string.Join("\n", testTransient.SimulatedValues)); 
		}
        [Test]
        public void TestCalculateTimeArray()
		{
            double[] timearray = Transient.CalculateTimeArray(0, 1000, 0.001);
            Console.WriteLine(string.Join("\n", timearray)); 
		}
        [Test]
        [TestCase(0.1)]
        [TestCase(1.0)]
        [TestCase(0.5)]
        public void TestClaculateIonDecayConstant(double value)
		{
            double[] timearray = Transient.CalculateTimeArray(0, 1000, 0.001);
            double result = Transient.CalculateIonDecayConstant(1000, value, 1000000);
            Console.WriteLine(result.ToString()); 
		}
        [Test]
        public void TestCalculateInstantaneousIonCount()
		{
            double[] timearray = Transient.CalculateTimeArray(0, 1, 0.001);
            int[] ionCountArray = Transient.CalculateInstantaneousIonCount(0.0, timearray, 100000);
            Console.WriteLine(string.Join("\n", ionCountArray)); 
		}
        [Test]
        public void TestCalculateInstantaneousTransient()
		{
            double[] result = Transient.SimulateTransient(3000, 2, 0.4, 9677, 0.01, 0.8, 0.000001, 10, 0.001, 1.0);
            Console.WriteLine(string.Join("\n", result)); 
		}
        [Test]
        public void TestCalculateInstantaneousTransientLargeProtein()
        {
            double[] result = Transient.SimulateTransient(100000, 40, 0.4, 9677, 0.01, 0.8, 0.000001, 10, 0.001, 1.0);
            Console.WriteLine(string.Join("\n", result));
        }

        [Test]
        public void TestSinFunction()
		{
            double[] timearray = Transient.CalculateTimeArray(0, 1, 0.001);
            double[] sinSeries = new double[timearray.Length];  
            for(int i = 0; i < timearray.Length; i++)
			{
                sinSeries[i] = Math.Sin(10 * timearray[i]);
            }
            Console.WriteLine(string.Join("\n", sinSeries)); 
		}
        [Test]
        public void TestGenerateIonCounts()
		{
            int[] testChargeStates = { 2, 3, 4, 5 };
            double[] testIsoInt = { 0.2, 0.3, 0.4, 0.1};
            double[] testIonCounts = { 2000,3000,4000,3000};
            double[] testMasses = { 2123, 2124, 2125, 2126};

            List<(int,double,double)> testResults = Transient.GenerateIonCounts(testChargeStates, testIsoInt, testIonCounts, testMasses); 
            foreach(var element in testResults)
			{
                Console.WriteLine(string.Join(" ;", element));
                Console.Write("\n"); 
			}

		}
        [Test]
        public void TestTransientGeneration()
		{
            int[] testChargeStates = { 2, 3, 4, 5 };
            double[] testIsoInt = { 0.2, 0.3, 0.4, 0.1 };
            double[] testIonCounts = { 2000, 3000, 4000, 3000 };
            double[] testMasses = { 2123, 2124, 2125, 2126 };

            List<(int, double, double)> testResults = Transient.GenerateIonCounts(testChargeStates, testIsoInt, testIonCounts, testMasses);
            InstrumentParameters intrPara = new(0.01, 0.95, 1.0, true, 0.65);
            List<Transient> transientResults = Transient.GenerateManyTransients(testResults, intrPara, 0.001, 10, 0.001, 0.6);
            double[] summedResults = Transient.SumMultipleTransients(transientResults);
            Console.WriteLine(string.Join("\n", transientResults[0].SimulatedValues)); 

        }
        [Test]
        public void TestTransientInitialiationII()
		{
            (int, double, double) testIonCountTuple = (2, 2123, 10000000);
            InstrumentParameters intrPara = new InstrumentParameters(0.01, 0.98, 300000, true, 0.95);
            Transient testTransient = new Transient(testIonCountTuple, intrPara, 0, 5, 0.001, false, 1.0);
            Console.WriteLine(string.Join("\n", testTransient.SimulatedValues)); 
		}
        [Test]
        public void CalculateFrequency()
		{
            double k = 300000;
            double q = 1 * 1.602E-19;
            double m = 55.934936 * 1.66054E-27;
            double result = Math.Sqrt((q/m) * k);
            Console.WriteLine(result.ToString());
		}
	}
}
