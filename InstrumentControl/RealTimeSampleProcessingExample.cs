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
using NUnit.Framework;
using UsefulProteomicsDatabases;
using Chemistry;
using System.IO;
using Proteomics;
using Proteomics.AminoAcidPolymer;
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
		public override Queue<MsDataScan> ScanProcessingQueue { get { return _scanProcessingQueue; } }
		public MS1DatabaseParser Database { get; }
		public MS1SearchEngine SearchEngine { get; }
		public int[] ScoreTable { get; private set; }
		public List<IsotopicEnvelope> Envelopes { get; private set; }
		private readonly int PeaksToKeep;
		public double[] BestMatchedPeaks { get; private set; } // May only be used for testing purpose, We shall see after the ICustomScan wrapper is built


		public RealTimeSampleProcessingExample(string databaseFileName, int ppmTolerance = 20, int peaksToKeep = 5)
        {
			Database = new(databaseFileName);
			SearchEngine = new(Database, new PpmTolerance(ppmTolerance));
			PeaksToKeep = peaksToKeep;
		}

		public RealTimeSampleProcessingExample(MS1DatabaseParser database, int ppmTolerance = 20, int peaksToKeep = 5)
        {
			Database = database;
			SearchEngine = new(Database, new PpmTolerance(ppmTolerance));
			PeaksToKeep = peaksToKeep;
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
			List<IsotopicEnvelope> envelopes = new();
            SearchEngine.PeakScorer(scan, out envelopes, out int[] scores);
            // TOTRY: passing in raw ScoreTable and Envelopes after wiping them above SearchEngine
            ScoreTable = scores;
			Envelopes = envelopes;



			// Takes the scored values unrecognizes them and adds them to a list of selected peaks, saving only the 5 peaks with the highest intensity
			// TODO Parallelize this
			// TOTRY: Adding all peaks to one list then pulling out the top 5 at the end, may have fewer operations and speed up processing overall
			List<IsotopicEnvelope> selectedPeaks = new();
			for (int i = 0; i < ScoreTable.Count(); i++)
            {
				if (ScoreTable[i] == 0)
                {
					if (selectedPeaks.Count() < PeaksToKeep)
                    {
						selectedPeaks.Add(Envelopes[i]);
                    }
					else if (selectedPeaks.Count() >= PeaksToKeep && Envelopes[i].TotalIntensity > selectedPeaks.Min(p => p.TotalIntensity))
                    {
						if (selectedPeaks.Any(p => SearchEngine.Tolerance.Within(Envelopes[i].MonoisotopicMass, p.MonoisotopicMass)))
						{
							continue;
                        }

						selectedPeaks.Remove(selectedPeaks.Where(p => p.TotalIntensity == selectedPeaks.Min(m => m.TotalIntensity)).Single());
						selectedPeaks.Add(Envelopes[i]);
                    }
                }
            }
			BestMatchedPeaks = selectedPeaks.Select(p => p.MonoisotopicMass).ToArray();
		}
	}
}

