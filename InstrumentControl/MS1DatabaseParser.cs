using Chemistry;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
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
	}
}
