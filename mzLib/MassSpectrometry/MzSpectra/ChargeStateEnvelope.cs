using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.MzSpectra
{
	public class ChargeStateEnvelope
	{
		public List<IsotopeCluster> IsotopeClusters { get; set; }
		public double MassGuess { get; set; }
		public double MonoisotopicMass { get; set; }
		public double AverageMass { get; set; }
		public double Score { get; set; }
		public double[] PossibleChargeStates { get; set; }
		public List<double> ChargeStateMatchesScores { get; set; }
		public double AdductMass { get; set; }

		public ChargeStateEnvelope(double adductMass = 1.007)
		{
			IsotopeClusters = new List<IsotopeCluster>();
			AdductMass = adductMass;
			MassGuess = 0;
		}
		public ChargeStateEnvelope(IsotopeCluster cluster)
		{
			IsotopeClusters = new List<IsotopeCluster>();
			AdductMass = cluster.AdductMass;
			MassGuess = 0;
			IsotopeClusters.Add(cluster);
		}
		public void CalculateMassGuessUsingIsotopomerSpacing()
		{
			double[] massGuessArray = new double[IsotopeClusters.Count];
			for (int i = 0; i < IsotopeClusters.Count; i++)
			{
				// cluster.ChargeState will have either a positive int charge or -1
				// if the charge state was unable to be calculated. 
				if (IsotopeClusters[i].ChargeState > 0)
				{
					massGuessArray[i] = IsotopeClusters[i].ChargeState * Math.Exp(IsotopeClusters[i].MaxPeakInCluster.Mz - Math.Log(AdductMass));
				}
			}
			MassGuess = massGuessArray.Average();
		}
		public void CaclulateMassGuessUsingChargeStateSpacing()
		{

		}
		public void UpdateMassGuess(double mz, int chargeState)
		{
			MassGuess += chargeState * Math.Exp(mz - Math.Log(AdductMass)) / (IsotopeClusters.Count + 1);
		}
		public void AddIsotopeCluster(IsotopeCluster cluster)
		{
			IsotopeClusters.Add(cluster);
			// update the mass guess if the charge state of the new cluster is valid and the mass guess has 
			// already been calculated initially
			if (MassGuess > 0 && cluster.ChargeState > 0)
			{
				UpdateMassGuess(cluster.MaxPeakInCluster.Mz, cluster.ChargeState);
			}
		}
		public void FindInitialIsotopicCluster(LogTransformedSpectra spectra, double peakFindingWindow)
		{
			// take top peak, grab all peaks within peakFindingWindowToleranceDa
			MzPeak topIntPeak = spectra.TransformedPeaksList
				.Aggregate((i1, i2) => i1.Intensity > i2.Intensity ? i1 : i2);

			var peaksInProposedCluster = FilterPeaks(spectra.TransformedPeaksList,
				topIntPeak.Mz, peakFindingWindow).ToList(); 
			// check for mass difference of an isotope and toss any that don't match. 
			// need to calculate the isotope mass difference still. 
			IsotopeCluster initialCluster = new();
			initialCluster.AddProposedPeaksAndFilterValidity(peaksInProposedCluster);
			IsotopeClusters.Add(initialCluster);
		}
		public void FindNextIsotopicCluster(IsotopeCluster cluster, LogTransformedSpectra spectra)
		{
			// next means z + 1, so the next cluster will be lower in m/z than the initial one
			// based on an intial IsotopicCluster, find all the other potential IsotopicClusters 
			// belonging to the same charge state envelope. 

			// start with the last computed isotope cluster.
			// if you don't get matched peaks with the last computed isotope cluster, then move to calculating above the first isotopecluster.

			if (cluster.ChargeState < 0)
			{
				// exits early if the charge state wasn't calculatable.
				return; 
			}
			// the next charge state is going to be +log(z1/z2)
			int nextChargeState = cluster.ChargeState + 1;
			double logMzOffset = Math.Log(cluster.ChargeState / nextChargeState);
			double nextMz = cluster.MaxPeakInCluster.Mz - logMzOffset;

			// check that the peak exists within a narrow range
			List<MzPeak> listOfMatchedPeaks = spectra.TransformedPeaksList.ToList().FindAll(i => PeakWithinRange(nextMz, i.Mz, 0.0001)); 
			if(listOfMatchedPeaks.Count == 0)
            {
				return; 
            }
			var filteredPeaks = FilterPeaks(spectra.TransformedPeaksList, nextMz, 0.05)
				.ToList();
			IsotopeCluster nextCluster = new();
			nextCluster.AddProposedPeaksAndFilterValidity(filteredPeaks);
			IsotopeClusters.Add(nextCluster); 
		}
		public bool PeakWithinRange(double initialPeak, double testPeak, double range)
		{
			double minVal = initialPeak + range;
			double maxVal = initialPeak - range;
			return testPeak < maxVal && testPeak > minVal;
		}
		public IEnumerable<MzPeak> FilterPeaks(List<MzPeak> peaksList, double mz, 
			double window)
        {; 
			double mzMin = mz - window;	
			double mzMax = mz + window;	
			return peaksList.Where(peak => peak.Mz >= mzMin && peak.Mz <= mzMax);	
        }

	}
}

