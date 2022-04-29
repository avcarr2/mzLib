using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.MzSpectra
{
	
	public class LogTransformedSpectra 
	{ 
		public List<MzPeak> TransformedPeaksList { get; set; }
		public double[] IntArray { get; set; }
		public double[] LogMzArray { get; set; }
		private double[] addcutMassLogMzSecondDerivative { get; set; }
		private double[,] chargeStateSeriesDiffs { get; set; }

		public LogTransformedSpectra(MzSpectrum spectrum, double adductMass)
		{
			TransformedPeaksList = spectrum.XArray
				.Zip(spectrum.YArray, (i, j) => new MzPeak(Math.Log(i) + Math.Log(adductMass), j))
				.ToList();
			IntArray = spectrum.XArray;
			LogMzArray = TransformedPeaksList.Select(i => i.Mz).ToArray(); 
		}
		
		public void UpdateSpectraOnFindingIsotopeCluster()
        {

        }
		
	}
	
}
