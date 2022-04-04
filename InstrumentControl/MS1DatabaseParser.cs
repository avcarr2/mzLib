using Chemistry;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
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
		public MS1DatabaseParser(string fileName)
		{
			// TODO: Add Ability to use other types of databases
			List<Protein> proteins = new();
			ProteinList = new();
			if (fileName.Contains(".fasta"))
			{
				proteins = ProteinDbLoader.LoadProteinFasta(fileName, true, DecoyType.None, false, out var dbErrors);
				DatabaseProcessing(proteins);
			}
			else
			{
				throw new Exception("Database file format not recognized");
			}
		}

		public MS1DatabaseParser(List<SimulatedProtein> proteinList)
		{
			ProteinList = proteinList.OrderBy(p => p.MonoisotopicMass).ToList();
			ProteinIndex = ProteinList.OrderBy(m => m.MonoisotopicMass).Select(p => p.MonoisotopicMass).ToArray();
		}

		public void DatabaseProcessing(List<Protein> proteins)
		{

			// TODO parallelize this
			// TODO add in ptm's (maybe gptmd style?)

			foreach (var protein in proteins)
			{
				SimulatedProtein simProt = new SimulatedProtein(protein);
				ProteinList.Add(simProt);
			}
			ProteinList = ProteinList.OrderBy(p => p.MonoisotopicMass).ToList();
			ProteinIndex = ProteinList.OrderBy(m => m.MonoisotopicMass).Select(p => p.MonoisotopicMass).ToArray();

		}

		// Takes in the filepath of a psmtsv and returns the proteins
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
	}
}
