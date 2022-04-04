using InstrumentControl;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
	public static class TestRealTimeSampleProcessingExample
	{
		[OneTimeSetUp]
		public static void Setup()
		{
			Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
		}

		[Test]
		// Test designed to get me into debug mode to evaluate what is going on within the code
		public static void TestRealTimeProcessing()
		{
			// Load in Database Data
			string fileName = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\Six_Protein_Standard.fasta";
			RealTimeSampleProcessingExample ProcessingExample = new(fileName, 20);


			// Load in msScans and send them through the scan processor
			//string filepath = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\2021-01-06_TopDownStandard_YeastFraction1.raw";
			string filepath = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\OneScanProteinStandard.mzML";

			List<MsDataScan> scans = new();
			string scanType = filepath.Split('.')[1].Trim();
			if (scanType.Equals("mzML"))
				scans = Mzml.LoadAllStaticData(filepath).GetAllScansList();
			if (scanType.Equals("raw"))
				scans = ThermoRawFileReader.LoadAllStaticData(filepath).GetAllScansList();
			scans = scans.Where(p => p.MsnOrder == 1).ToList();

			foreach (var scan in scans)
			{
				
				ProcessingExample.ScanProcessingQueue.Enqueue(scan);
				ProcessingExample.ProteoformProcessingEngine();

			}

			Assert.That(true);
		}

		/// <summary>
		/// Test takes a MetaMorpheus TopDown search and removes a third of the protein results then runs it thorugh the RealTimeSampleProcessingExample with the non-removed proteins as the database
		/// Returns true if at least 70% of the removed database proteins are selected for fragmentation 
		/// </summary>
		[Test]
		public static void TestingRealTimeProcessingOnSearchResults()
        {
			double percentToRemove = 30;

			// Loads in scans
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			//string filepath = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\2021-01-06_TopDownStandard_YeastFraction1.raw"; // for me to run with full LCMS run
			var scans = Mzml.LoadAllStaticData(filepath).GetAllScansList();

			// Loads in MM-TD search results of the above scans, pulls out the top scoring 100, and treats half as the database
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			//string psmFile = @"C:\Users\Nic\Desktop\FileAccessFolder\API\RealTimeProcessingExample\TDYestFractionMMFullResults\Task1-SearchTask\AllProteoforms.psmtsv"; //for me to run with full LCMS results
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			List<SimulatedProtein> removedProteins = proteins.GetRange(0, (int)(proteins.Count() * (percentToRemove / 100)));
			proteins.RemoveRange(0, (int)(proteins.Count() * (percentToRemove / 100)));

			MS1DatabaseParser database = new MS1DatabaseParser(proteins);
			RealTimeSampleProcessingExample processingExample = new(database);

			// Send the scans through the search engine and get a list of what it chooses to fragment
			List<double> selectedMasses = new();
			foreach (var scan in scans)
			{
				processingExample.ScanProcessingQueue.Enqueue(scan);
				processingExample.ProteoformProcessingEngine();
				selectedMasses.AddRange(processingExample.BestMatchedPeaks);
			}
			var groupedmasses = selectedMasses.GroupBy(p => Math.Round(p, 2)).ToList();

			// Iterates through all proteins that were removed from the database, creates an array of their masses, and sets the mass to 0 if that mass was identified to be fragmented by the program
			double[] removedProteinMasses = removedProteins.Select(p => p.MonoisotopicMass).ToArray();
            for (int i = 0; i < removedProteins.Count(); i++)
            {
				if (selectedMasses.Any(p => processingExample.SearchEngine.Tolerance.Within(p, removedProteins[i].MonoisotopicMass)))
				{
					removedProteinMasses[i] = 0;
				}
			}

			double percentOfRemovedProteins = (double)removedProteins.Count() * 0.30;
			int countOfRemovedProteinsThatWereNotSelected = removedProteinMasses.Count(p => p != 0);
			Assert.GreaterOrEqual(percentOfRemovedProteins, countOfRemovedProteinsThatWereNotSelected);

		}

		[Test]
		public static void TestTimingOfProteoformProcessingEngine()
        {
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			var scans = Mzml.LoadAllStaticData(filepath).GetAllScansList();
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			MS1DatabaseParser database = new MS1DatabaseParser(proteins);
			RealTimeSampleProcessingExample processingExample = new(database);

			var stopWatch = new Stopwatch();
			List<double> times = new();
			foreach (var scan in scans)
            {
				stopWatch.Start();
				processingExample.ScanProcessingQueue.Enqueue(scan);
				processingExample.ProteoformProcessingEngine();
				stopWatch.Stop();
				times.Add(stopWatch.Elapsed.TotalMilliseconds);
				stopWatch.Reset();
			}
			double minTime = times.Min();
			double maxTime = times.Max();
			double avgTime = times.Average();

			Console.WriteLine("Average Time: {0} /n Max Time: {1} /n Min Time: {2}", avgTime, maxTime, minTime);
			Assert.That(250 >= avgTime);
			Assert.That(250 >= maxTime);
		}
	}

}
