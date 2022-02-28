using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InstrumentControl
{
	public static class MSScanContainerNames
	{
		public static readonly HashSet<string> MSScanContainerKeys; 
		static MSScanContainerNames()
		{
			// may need to add more keys in the future. 
			MSScanContainerKeys = new HashSet<string>()
			{
				"Total Ion Current",
				"Scan Low Mass",
				"Scan High Mass",
				"Scan Start time (min)",
				"Scan Number",
				"Base Peak Intensity",
				"Base Peak Mass",
				"Scan Mode",
				"AGC",
				"Micro Scan Count",
				"Ion Injection Time (ms)",
				"Elapsed Scan Time (sec)",
				"Average Scan by Inst",
				"Orbitrap Resolution",
				"API Process Delay",
				"Dependency Type",
				"Access Type",
				"Access ID",
				"Conversion Parameter I",
				"Conversion Parameter A",
				"Conversion Parameter B",
				"Conversion Parameter C",
				"Conversion Parameter D",
				"Conversion Parameter E",
				"Temperature Comp. (ppm)",
				"RF Comp. (ppm)",
				"Space Charge Comp. (ppm)",
				"Resolution Comp. (ppm)",
				"Number of LM Found",
				"LM Correction (ppm)",
				"RawOvFtT",
				"Injection t0",
				"Reagent Ion Injection Time (ms)",
				"FAIMS Voltage On",
				"FAIMS CV",
				"Multi Inject Info",
				"Master Scan Number",
				"Monoisotopic M/Z",
				"Charge State",
				"Error in isotopic envelope fit",
				"HCD Energy",
				"HCD Energy eV",
				"MS2 Isolation Width",
				"SPS Masses",
				"SPS Masses Continued"
			}; 
		}
	}
}
