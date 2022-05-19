using InstrumentControl;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
	public static class TestInstrumentControl
	{
		[OneTimeSetUp]
		public static void Setup()
		{
			Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
		}

		/// <summary>
		/// Test takes a MetaMorpheus TopDown search and removes a third of the protein results then runs it thorugh the RealTimeSampleProcessingExample with the non-removed proteins as the database
		/// Returns true if at least 70% of the removed database proteins are selected for fragmentation 
		/// </summary>
		[Test]
		public static void TestingRealTimeProcessingOnSearchResults()
        {
			double percentToRemove = 30;
			double percentToMatch = 80;
			int threads = 8;

			// Loads in scans
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			List<MsDataScan> scans = MS1DatabaseParser.LoadAllScansFromFile(filepath);

			// Loads in MM-TD search results of the above scans, pulls out the top scoring 100, and treats half as the database
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			List<SimulatedProtein> removedProteins = proteins.GetRange(0, (int)(proteins.Count() * (percentToRemove / 100)));
			proteins.RemoveRange(0, (int)(proteins.Count() * (percentToRemove / 100)));

			MS1DatabaseParser database = new MS1DatabaseParser(proteins);
			RealTimeSampleProcessingExample processingExample = new(database, threads);

			// Send the scans through the search engine and get a list of what it chooses to fragment
			List<double> selectedMasses = new();
			foreach (var scan in scans)
			{
				processingExample.ScanProcessingQueue.Enqueue(scan);
				processingExample.ProteoformProcessingEngine();
				selectedMasses.AddRange(processingExample.PeaksToFragment);
			}

			// Iterates through all proteins that were removed from the database, creates an array of their masses, and sets the mass to 0 if that mass was identified to be fragmented by the program
			double[] removedProteinMasses = removedProteins.Select(p => p.MonoisotopicMass).ToArray();
            for (int i = 0; i < removedProteins.Count(); i++)
            {
				if (selectedMasses.Any(p => processingExample.SearchEngine.Tolerance.Within(p, removedProteins[i].MonoisotopicMass)))
				{
					removedProteinMasses[i] = 0;
				}
			}

			double percentOfRemovedProteinsSelectedForFragmentation = removedProteinMasses.Count(p => p == 0) / (double)removedProteins.Count * 100;
			Assert.GreaterOrEqual(percentOfRemovedProteinsSelectedForFragmentation, percentToMatch);
		}

		[Test]
		public static void TestTimingOfProteoformProcessingEngine()
        {
			int threads = 20;
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			List<MsDataScan> scans = MS1DatabaseParser.LoadAllScansFromFile(filepath);
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			MS1DatabaseParser database = new MS1DatabaseParser(proteins, threads);
			RealTimeSampleProcessingExample processingExample = new(database, threads);

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

			Console.WriteLine("Average Time: {0} \n Max Time: {1} \n Min Time: {2}", times.Average(), times.Max(), times.Min());
			Assert.That(250 >= times.Average());
			Assert.That(250 >= times.Max());
		}

		[Test]
		public static void TestMS1SearchEngineFindPeakWithinDatabase()
        {
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			List<MsDataScan> scans = MS1DatabaseParser.LoadSelectScansFromFile(filepath, 32, 34);
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			MS1DatabaseParser database = new MS1DatabaseParser(proteins);
			RealTimeSampleProcessingExample processingExample = new(database, 6);

			for (int i = 0; i < scans.Count; i++)
            {
				processingExample.ScanProcessingQueue.Enqueue(scans[i]);
				processingExample.ProteoformProcessingEngine();
				short[] testTable = new short[processingExample.SearchEngine.ScoreTable.Length];
				if (i == 0)
					testTable[42] = 1;
				else if (i == 1)
					testTable[51] = 1;
				Assert.AreEqual(testTable, processingExample.SearchEngine.ScoreTable);
			}
		}

		[Test]
		public static void TestMS1SearchEngineMultiThreading()
        {
			string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMS1.mzML");
			List<MsDataScan> scans = MS1DatabaseParser.LoadSelectScansFromFile(filepath, 32);
			string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\TDYeastFractionMMResult.psmtsv");
			List<SimulatedProtein> proteins = MS1DatabaseParser.GetSimulatedProteinsFromPsmtsv(psmFile);
			MS1DatabaseParser database = new MS1DatabaseParser(proteins);
			RealTimeSampleProcessingExample processingExample = new(database, 1, 20, 5);
			RealTimeSampleProcessingExample processingExample2 = new(database, 2, 20, 5);
			RealTimeSampleProcessingExample processingExample6 = new(database, 6, 20, 5);

			processingExample.ScanProcessingQueue.Enqueue(scans[0]);
			processingExample2.ScanProcessingQueue.Enqueue(scans[0]);
			processingExample6.ScanProcessingQueue.Enqueue(scans[0]);
			processingExample.ProteoformProcessingEngine();
			processingExample2.ProteoformProcessingEngine();
			processingExample6.ProteoformProcessingEngine();

			Assert.AreEqual(processingExample.PeaksToFragment, processingExample2.PeaksToFragment);
			Assert.AreEqual(processingExample.PeaksToFragment, processingExample6.PeaksToFragment);
			Assert.AreEqual(processingExample.SearchEngine.ScoreTable, processingExample2.SearchEngine.ScoreTable);
			Assert.AreEqual(processingExample.SearchEngine.ScoreTable, processingExample6.SearchEngine.ScoreTable);
		}

	}

}
