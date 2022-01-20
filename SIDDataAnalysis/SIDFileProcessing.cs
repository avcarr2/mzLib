using IO.MzML;
using MassSpectrometry;
using System;
using System.Diagnostics;
using System.IO;
using IO.ThermoRawFileReader;
using System.Linq;
using MzLibUtil;
using System.Collections.Generic;

namespace SIDDataAnalysis
{
	public class SIDFileProcessing
	{
		// should've probably just made these methods in MSDataScan,but whatever
		// changes order of sidMS scan to MS2 in a MS1/sidMS2 experimental format 
		public static List<MsDataScan> ConvertMS1Heading(List<MsDataScan> scansList)
		{
			for (int i = 0; i < scansList.Count; i++)
			{
				if (scansList[i].OneBasedScanNumber % 2 == 1)
				{
					continue;
				}
				else
				{
					int precursorScanNumber = scansList[i].OneBasedScanNumber - 1;
					// adds necessary MS2 information. Without, file writing throws an error. 
					scansList[i].SetMsnOrder(2);
					scansList[i].SetOneBasedPrecursorScanNumber(precursorScanNumber);
					scansList[i].SetIsolationMz(1000);
					scansList[i].SetIsolationWidth(2000);
					scansList[i].SetSelectedIonMz(1000);
				}
			}
			return scansList;
		}
		public static MzSpectrum DeconvoluteSpectrum(MsDataScan scan)
		{
			// create a new MzRange to pass to .Deconvolute(); 
			MzRange range = new(scan.MassSpectrum.XArray.Min(), scan.MassSpectrum.XArray.Max());
			// constant parameters to pass to .Deconvolute(); 
			int minCharge = 1;
			int maxCharge = 60;
			double deconTolerance = 4.00;
			double intensityRatioLimit = 100;
			List<IsotopicEnvelope> deconIE = scan.MassSpectrum.Deconvolute(range, minCharge, maxCharge, deconTolerance, intensityRatioLimit).ToList();
			// I think dictionary is easiest here, but I don't think it's the fastest solution. 
			Dictionary<double, double> monoMassAndIntDict = new Dictionary<double, double>();
			double[] uniqueMono = deconIE.Select(i => i.MonoisotopicMass).Distinct().ToArray();
			// For each unique monoIsotopic mass, this loop finds all matching uniqueMono[i] and 
			// sums the intensity associated with that value. 
			for (int i = 0; i < uniqueMono.Length; i++)
			{
				monoMassAndIntDict.Add(uniqueMono[i], deconIE.FindAll(delegate (IsotopicEnvelope ie)
				{
					return ie.MonoisotopicMass == uniqueMono[i]; // delegate function to return all monoisotopic 
					// masses in List<IsotopicEnvelope> that are equal to the current uniqueMono.  
				}).Select(i => i.TotalIntensity) // select the total intensity
				.Sum()); // sum the total intensities. 
			}
			// return a new MzSpectrum that contains the new peaks. 
			return  new(monoMassAndIntDict.Keys.ToArray(), monoMassAndIntDict.Values.ToArray(), true);
		}
		public static List<MsDataScan> DeconvoluteAllSpectra(List<MsDataScan> scans) 
		{ 
			foreach(MsDataScan scan in scans)
			{
				scan.OverwriteMzSpectrum(DeconvoluteSpectrum(scan)); 
			}
			return scans; 
		} 
		public static List<MsDataScan> LoadAndReadSIDMSExperiment(string filePath)
		{
			return ThermoRawFileReader.LoadAllStaticData(filePath, null, 1).GetAllScansList(); 
		}
		public static List<MsDataScan> DeconvoluteAndUpdateHeaders(List<MsDataScan> scans)
		{
			// Must deconvolute while all scans are still MS1s! Deconvolution is different for MS2s!
			return ConvertMS1Heading(DeconvoluteAllSpectra(scans)); 
		}
	}
}