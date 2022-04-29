using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.MzSpectra
{
	public class IsotopeCluster
	{
		public List<MzPeak> PeaksInClusterList { get; set; }
		public int ChargeState { get; set; }
		public MzPeak MaxPeakInCluster { get; set; }
		public int NumberPeaksInCluster { get; set; }
		public double MonoisotopicMass { get; set; }
		public double AverageMass { get; set; }
		public double Score { get; set; }
		public double CumulativeIntensity { get; set; }
		public double AdductMass { get; set; }
		
		// difference between an isotope: 
		// 1/EXP(\delta(LOG(MZ)) - LOG(adduct mass)) = isotope mass.
		// if the calculation equals the isotope mass, the peaks are isotopomers. 
		// this will require centroided data, because profile data has too many peaks
		// and the code will run slowly. 
		public IsotopeCluster(double adductMass = 1.007)
		{
			PeaksInClusterList = new List<MzPeak>();
			MaxPeakInCluster = new MzPeak(0, 0);
			AdductMass = adductMass;
		}
		public IsotopeCluster(List<MzPeak> listPeaks, double adductMass = 1.007)
		{
			PeaksInClusterList = listPeaks;
			MaxPeakInCluster = listPeaks.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2);
			AdductMass = adductMass;
			NumberPeaksInCluster = PeaksInClusterList.Count;
			CalculateChargeState(); 
		}
		public void AddProposedPeaksAndFilterValidity(List<MzPeak> proposedCluster)
        {
			// This method currently works, however, it's only as good as your peak centroiding.
			// And it's also greedy, i.e. it will select all peaks that yield the isotope mass to 
			// be a part of the IsotopeCluster.
			// get all pairs -> compute that a peak is an isotope mass away 
			var validPeakPairsList = GetAllPairs(proposedCluster)
				.ToList()
				.FindAll(i => PeakIsIsotopeMassAway(i.Item1.Mz, i.Item2.Mz, isotopeMass: 1.003,
				tolerance: 0.0002, adductMass: 1.007)); // output is a list of mzpeakpairs.

			// create the final list of unique MzPeaks that are part of the isotopic envelope. 
			List<MzPeak> flattenedList = new();
			foreach (var pair in validPeakPairsList)
			{
				flattenedList.Add(pair.Item1);
				flattenedList.Add(pair.Item2);
			}



			PeaksInClusterList = flattenedList.Distinct().ToList();
			MaxPeakInCluster = flattenedList.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2);
			NumberPeaksInCluster = PeaksInClusterList.Count;
			CalculateChargeState(); 
		}
		public void RemoveExtraPeaksFromIsotopicEnvelope(List<MzPeak> peaks)
        {

        }
		public static IEnumerable<(T, T)> GetAllPairs<T>(IList<T> source)
		{
			return source.SelectMany((_, i) => source.Where((_, j) => i < j),
				(x, y) => (x, y));
		}
		public static bool PeakIsIsotopeMassAway(double peak1, double peak2, 
			double isotopeMass, double tolerance, double adductMass)
		{
			double diff = peak2 - peak1;
			// create a switch to handle cases
			switch (diff)
			{
				case var _ when diff.Equals(0): return false;
				case var _ when diff < 0: _ = Math.Abs(diff); break;
				default: break;
			}
			double isotopeMassTest = IsotopomerCalculation(diff,isotopeMass, adductMass);
			double maxVal = isotopeMass + tolerance; 
			double minVal = isotopeMass - tolerance;
			return isotopeMassTest >= minVal && isotopeMassTest <= maxVal;
		}
		public static double IsotopomerCalculation(double logMzDiff,double isotopeMass, double adductMass)
		{
			double lnAdductMass = Math.Log(adductMass);
			double lnIsotopeMass = Math.Log(isotopeMass);
			return 1 / Math.Exp(logMzDiff - lnAdductMass + lnIsotopeMass);
		}
		public bool CheckValidPeak(MzPeak peak)
		{
			return true;
		}
		public bool AddIsotopicPeak(MzPeak peak)
		{
			if (!CheckValidPeak(peak))
			{
				return false;
			}
			PeaksInClusterList.Add(peak);
			// update the cumulative intensity
			CumulativeIntensity += peak.Intensity;
			NumberPeaksInCluster += 1;
			CheckForNewMaxIntensity(peak);
			UpdateScore();
			return true;
		}
		public double UpdateScore()
		{
			return 0.0;
		}
		public void CalculateAverageMass() { }
		public void CalculateMonoisotopicMass() { }
		public void CheckForNewMaxIntensity(MzPeak peak)
		{
			if (peak.Intensity > MaxPeakInCluster.Intensity)
			{
				MaxPeakInCluster = peak;
			}
		}
		public void CalculateChargeState()
		{
			// order by m/z.
			double[] mzVals = PeaksInClusterList
				.OrderBy(i => i.Mz)
				.Select(i => i.Mz)
				.ToArray();

			// charge of the ion = adductMass / (exp(mzVals[i]) - exp(mzVals[i + 1])
			double[] charge = new double[mzVals.Length - 1];
			for (int i = 0; i < charge.Length; i++)
			{
				charge[i] = AdductMass / (Math.Exp(mzVals[i + 1]) - Math.Exp(mzVals[i]));
			}
			CheckValidChargeState(charge, 0.01, out int intCharge);
			ChargeState = intCharge;
		}
		public static bool CheckValidChargeState(double[] chargeStateGuesses, double stdevThresh,
			out int intCharge)
		{
			// returns charge state as integer value if the stdev of the charge state guess is below 
			// a threshold. 
			// returns -1 otherwise. 
			// The -1 might mean that there's bad isotopes in the isotopeCluster, 
			// or the isotopeCluster is not isotopically resolved. 
			double mean = chargeStateGuesses.Average();
			double stdev = Math.Sqrt(chargeStateGuesses.Average(v => Math.Pow(v - mean, 2)));
			if (stdev <= stdevThresh)
			{
				intCharge = (int)Math.Round(mean);
				return true;
			}
			intCharge = -1;
			return false;
		}
	}
	
}

