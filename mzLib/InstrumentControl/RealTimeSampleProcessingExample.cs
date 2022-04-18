using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Thermo.Interfaces.FusionAccess_V1;
using Thermo.Interfaces.FusionAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.Control.Scans;
using Thermo.Interfaces.SpectrumFormat_V1;
using Thermo.Interfaces.FusionAccess_V1.Control.Scans;
using System.Threading;
using MassSpectrometry;
using MzLibUtil;

namespace InstrumentControl
{
	public class RealTimeSampleProcessingExample : DataReceiver
	{
		/*
		 Method classes inheriting from DataReceiver must implement a ScanProcessingQueue
		(Queue<MsDataScan>) and the MSScanContainer_MsScanArrived method.

		To implement the property and method, you must use the override keyword.

		After the abstract method overrides, implement the "engine" method, where 

		 */

		// private field used to initialize the scan processing queue
		private Queue<MsDataScan> _scanProcessingQueue = new Queue<MsDataScan>();
		// ScanProcessingQueue is the public property; it uses the getter to set the private
		// field to the public property. This is the correct way to initialize a publicly
		// available object like Queues, Lists, and Dicts as properties.
		private readonly int PeaksToKeep;
		public override Queue<MsDataScan> ScanProcessingQueue { get { return _scanProcessingQueue; } }
		public MS1DatabaseParser Database { get; }
		public MS1SearchEngine SearchEngine { get; }
		public List<IsotopicEnvelope> Envelopes { get; private set; }
		public double[] PeaksToFragment { get; private set; } 


		public RealTimeSampleProcessingExample(string databaseFileName, int threads, int ppmTolerance = 20, int peaksToKeep = 5)
        {
			Database = new(databaseFileName, threads);
			SearchEngine = new(Database, new PpmTolerance(ppmTolerance), threads);
			PeaksToKeep = peaksToKeep;
			PeaksToFragment = new double[peaksToKeep];
		}

		public RealTimeSampleProcessingExample(MS1DatabaseParser database, int threads, int ppmTolerance = 20, int peaksToKeep = 5)
        {
			Database = database;
			SearchEngine = new(Database, new PpmTolerance(ppmTolerance), threads);
			PeaksToKeep = peaksToKeep;
			PeaksToFragment = new double[peaksToKeep];
		}

		public RealTimeSampleProcessingExample() { }

		public override void MSScanContainer_MsScanArrived(object sender, MsScanEventArgs e)
		{
			// put scan in queue
			// Keeping IMScan open will hog the memory resources, so you need to get the  

			using (IMsScan scan = e.GetScan())
			{
				// need to quickly convert to something else
				// ScanProcessingQueue.Enqueue(scan);
				ScanProcessingQueue.Enqueue(new MsDataScan(scan));
			}
			// Perform data processing below here. For clarity, use the minimum
			// number of method calls possible, preferablly a single method,
			// for example, ProteoformIdentifcationEngine(ProteoformIDParameters params);
			// Engine should return either void or an ICustomScan object, but I need to implement the
			// code to facilitate the building and sending of the ICustomScan. 

			ProteoformProcessingEngine();

			// TODO: Implement ICustomScan wrapper to facilitate custom scan sending.
			// May need to wait till Chris Rose has figured everything with Thermo out. 
		}

		// other methods used in processing go below this comment:
		public void ProteoformProcessingEngine()
		{
			// engine should try to return void, bool, or an ICustomScan object (or both bool and ICustomScan) 

			// Deque the scan to be analyzed and clear any data from previous scans entereing this method
			MsDataScan scan = ScanProcessingQueue.Dequeue();

            List<IsotopicEnvelope> envelopes = SearchEngine.FindPeakWithinDatabase(scan);
			// Could assign scores to ScoreTable and envelopes to Envelopes but I dont see the utility of them being accessible outside of the method atm
            List<IsotopicEnvelope> peaksNotFoundInDatabase = envelopes.Where(n => SearchEngine.ScoreTable[envelopes.IndexOf(n)] == 0).OrderByDescending(p => p.TotalIntensity).ToList();
            PeaksToFragment = peaksNotFoundInDatabase.Take(PeaksToKeep).Select(p => p.MonoisotopicMass).ToArray();
		}
	}
}

