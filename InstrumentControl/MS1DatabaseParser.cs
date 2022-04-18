using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace InstrumentControl

{	/// <summary>
	/// Class designed to process to process databases
	/// </summary>
    public class MS1DatabaseParser
    {
		public List<SimulatedProtein> ProteinList { get; private set; }
		public double[] ProteinIndex { get; private set; }
		private int Threads;

		/// <summary>
		/// Construct for MS1DatabaseParser that loads in the database from its file location.
		/// Currently only supports fasta files
		/// </summary>
		public MS1DatabaseParser(string fileName, int threads = 1)
		{
			// TODO: Add Ability to use other types of databases
			List<Protein> proteins = new();
			ProteinList = new();
			Threads = threads;
			if (fileName.Contains(".fasta"))
			{
				proteins = ProteinDbLoader.LoadProteinFasta(fileName, true, DecoyType.None, false, out var dbErrors);
				LoadProteinDataAsSimulatedProtein(proteins);
			}
			else
			{
				throw new Exception("Database file format not recognized");
			}
		}

		/// <summary>
		/// Constructor for MS1DatabaseParser that accepts a list of SimulatedProteins
		/// </summary>
		/// <param name="proteinList"></param>
		public MS1DatabaseParser(List<SimulatedProtein> proteinList, int threads = 1)
		{
			Threads = threads;
			ProteinList = proteinList.OrderBy(p => p.MonoisotopicMass).ToList();
			ProteinIndex = ProteinList.OrderBy(m => m.MonoisotopicMass).Select(p => p.MonoisotopicMass).ToArray();
		}

		/// <summary>
		/// Loads in proteins as SimulatedProteins and populates the respective method fields
		/// </summary>
		public void LoadProteinDataAsSimulatedProtein(List<Protein> proteins)
		{
			// TODO: Add in PTM's
			int[] threads = Enumerable.Range(0, Threads).ToArray();
			Parallel.ForEach(threads, index =>
			{
				for (; index < proteins.Count; index += Threads)
				{
					SimulatedProtein simProt = new(proteins[index]);
					lock (ProteinList)
                    {
						ProteinList.Add(simProt);
					}
				}
			});

			//foreach (var protein in proteins)
			//{
			//	SimulatedProtein simProt = new SimulatedProtein(protein);
			//	ProteinList.Add(simProt);
			//}
			ProteinList = ProteinList.OrderBy(p => p.MonoisotopicMass).ToList();
			ProteinIndex = ProteinList.OrderBy(m => m.MonoisotopicMass).Select(p => p.MonoisotopicMass).ToArray();
		}

		/// <summary>
		/// Takes in the filepath of a psmtsv and returns them as a list of SimulatedProteins.
		/// Currently only returns those who do not have any ptm's
		/// </summary>
		public static List<SimulatedProtein> GetSimulatedProteinsFromPsmtsv(string filepath, bool onlyBaseSequences = true)
        {
			List<SimulatedProtein> proteins = new();
			foreach (string line in File.ReadAllLines(filepath))
			{
				var split = line.Split(new char[] { '\t' });
				if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
				{
					continue;
				}
				string baseSequence = split[12];
				string fullSequence = split[13];
				string accession = split[25];

				// If there are no PTM's as this is not implemented yet
				if (onlyBaseSequences && fullSequence.Equals(baseSequence))
				{
					SimulatedProtein protein = new(new Protein(baseSequence, accession));
					proteins.Add(protein);
				}
				// If you account for PTM's
				if (!onlyBaseSequences)
                {
					SimulatedProtein protein = new(new Protein(baseSequence, accession));
					proteins.Add(protein);
                }
			}
			return proteins;
		}

		/// <summary>
		/// Creates a List of MsDataScans from a spectra file. Currently supports MzML and raw
		/// </summary>
		public static List<MsDataScan> LoadAllScansFromFile(string filepath)
        {
			List<MsDataScan> scans = new();
			if (filepath.EndsWith(".mzML"))
				scans = Mzml.LoadAllStaticData(filepath).GetAllScansList();
			else if (filepath.EndsWith(".raw"))
				scans = ThermoRawFileReader.LoadAllStaticData(filepath).GetAllScansList();
			else
            {
				throw new MzLibException("Cannot load spectra");
            }
			return scans;
		}

		/// <summary>
		/// Creates a List of MsDataScans from a spectra file. Takes the OneBasedScanNumbers as inputs. 
		/// </summary>
		public static List<MsDataScan> LoadSelectScansFromFile(string filepath, int start, int end = -1)
        {
			if (end == -1)
            {
				end = start + 1;
            }
			List<MsDataScan> scans = LoadAllScansFromFile(filepath);
			List<MsDataScan> trimmedScans = new();
			//        for (; start < end + 1; start++)
			//        {
			//trimmedScans.Add(scans[start - 1]);
			//        }
			trimmedScans = scans.GetRange(start - 1, (end - start));
			return trimmedScans;

		}

	}
}
