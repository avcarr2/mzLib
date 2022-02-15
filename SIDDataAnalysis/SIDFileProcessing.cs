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
			int minCharge = 12;
			int maxCharge = 80;
			double deconTolerance = 4.00;
			double intensityRatioLimit = 1000000000;
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
			List<MsDataScan> scans = ThermoRawFileReader.LoadAllStaticData(filePath).GetAllScansList();
			return scans; 
		}
		public static List<MsDataScan> DeconvoluteAndUpdateHeaders(List<MsDataScan> scans)
		{
			// Must deconvolute while all scans are still MS1s! Deconvolution is different for MS2s!
			return ConvertMS1Heading(DeconvoluteAllSpectra(scans)); 
		}

		// Need to write the code to convert every three scans into a single scan. 
		// May need to normalize it? Use tic to normalize. 
		public static MsDataScan SumSpectra(MsDataScan parentScan, List<MsDataScan> daughterScansToSum, double minMZ, double maxMZ)
		{
			// ** will need to reorder the entire spectra. **
			// create new scan that will be formed from the old scans
			// normalize the scans according to TIC
			// add each sub scan 
			// return the subscan

			// NOTE: for this part, assign the one based scan number and the precursor scan number to the parentScan, 
			// then assign the new one based for the entire scan in a separate method

			// Get each mz and intensities values for each scan
			List<MzPeak> mzPeaks = new();
			// Need to calculate the new total ion current for the composited scan
			double cumulativeTIC = 0; 
			foreach(MsDataScan scan in daughterScansToSum)
			{
				// extract all intensity vs mz values into a list of MzPeaks. 
				// MzPeaks contains property mz and intensity that make it more 
				// amenable to combining the spectra. This method avoids nested for loops 
				// and relying on both double[,] and dictionaries and their interconversion
				double scanMin = scan.MassSpectrum.XArray.Min();
				double scanMax = scan.MassSpectrum.XArray.Max();

				cumulativeTIC += scan.TotalIonCurrent; 
				mzPeaks.AddRange(scan.MassSpectrum.Extract(scanMin, scanMax).ToList());
			}

			// need to calculate the selected ion and isolation width
			double compositeScanSelectedIonMz = (maxMZ - minMZ) / 2;
			double compositeScanIsolationWidthMZ = maxMZ - minMZ; 
			// break the List<MzPeaks> into two arrays for passing to MzSpectrum creation. 
			double[] xArray = mzPeaks.Select(i => i.Mz).ToArray();
			double[] yArray = mzPeaks.Select(i => i.Intensity).ToArray(); 
			
			// create new composite spectrum and fill it with mz and intensity values

			// Required information not obtainable from the parent scan properties
			MzSpectrum compositeSpectrum = new MzSpectrum(xArray, yArray, false);
			int compositeScanOrder = 2;
			MzRange compositeScanRange = new(compositeSpectrum.XArray.Min(), compositeSpectrum.XArray.Max());
			string compositeScanFilter = @"FTMS + p NSI sid=60.00 Full ms [500.0000-2000.0000]";
			// paceholder NativeID to be replaced later in the workflow by a different method. 
			string compositeScanNativeID = @"controllerType=0\tcontrollerNumber=1\tscan=" + (parentScan.OneBasedPrecursorScanNumber + 1).ToString(); 

			// create new scan: 
			MsDataScan compositeScan = new MsDataScan(compositeSpectrum,
				oneBasedScanNumber: parentScan.OneBasedScanNumber + 1,
				msnOrder: compositeScanOrder,
				isCentroid: false,
				polarity: parentScan.Polarity,
				// use first daughter scan retention time
				retentionTime: daughterScansToSum[0].RetentionTime,
				scanWindowRange: compositeScanRange,
				scanFilter: compositeScanFilter,
				// mzAnalyzer may be different in future experiments. 
				mzAnalyzer: parentScan.MzAnalyzer,
				totalIonCurrent: cumulativeTIC,
				// calculable if in the future, we want unbroken TICs viewable in seeMS or similar program. 
				injectionTime: null,
				noiseData: null,
				nativeId: compositeScanNativeID,
				selectedIonMz: compositeScanSelectedIonMz,
				selectedIonChargeStateGuess: null,
				selectedIonIntensity: null,
				isolationMZ: compositeScanSelectedIonMz,
				isolationWidth: compositeScanIsolationWidthMZ,
				dissociationType: DissociationType.Unknown,
				oneBasedPrecursorScanNumber: parentScan.OneBasedPrecursorScanNumber,
				selectedIonMonoisotopicGuessMz: null, 
				hcdEnergy: null);

			return (compositeScan); 
		}
		public static List<MsDataScan> SumSpectraAcrossAllScans(List<MsDataScan> allScans, int numberDaughterScans)
		{
			List<MsDataScan> processedScans = new(); 
			int i = 0; 
			while(i < allScans.Count - numberDaughterScans)
			{
				MsDataScan parentScan = allScans[i];
				List<MsDataScan> daughterScans = new();
				for (int j = i + 1; j < (i + numberDaughterScans + 1); j++)
				{
					daughterScans.Add(allScans[j]); 
				}
				double minMZ = daughterScans[0].MassSpectrum.XArray.Min();
				double maxMZ = daughterScans[numberDaughterScans - 1].MassSpectrum.XArray.Max();

				MsDataScan compositedScan = SumSpectra(parentScan, daughterScans, minMZ, maxMZ);

				processedScans.Add(parentScan);
				processedScans.Add(compositedScan); 

				daughterScans.Clear(); 
				i += numberDaughterScans + 1; 
			}

			List<MsDataScan> updatedScanNumbersScanList = UpdateScanAndPrecursorNumber(processedScans); 

			return updatedScanNumbersScanList; 
		}

		// Definitely need to write a test for this
		public static List<MsDataScan> UpdateScanAndPrecursorNumber(List<MsDataScan> processedScans)
		{
			// change all the one-based spectrum numbers

			// if scan MsnOrder = 1, need to give it the new one-based scan number and update the precursor scan number
			// if scan MsnOrder = 2, need to give it the new one-based scan number and 
			int scanNumber = 1;
			foreach (MsDataScan scan in processedScans)
			{
				if (scan.MsnOrder == 1)
				{
					scan.SetOneBasedScanNumber(scanNumber);
					scan.SetNativeID("controllerType=0 controllerNumber=1" + " scan=" + scanNumber.ToString()); 
					scanNumber++;
					continue;
				}
				else
				{
					scan.SetOneBasedScanNumber(scanNumber); 
					scan.SetOneBasedPrecursorScanNumber(scanNumber - 1);
					scan.SetNativeID("controllerType=0 controllerNumber=1" + " scan=" + scanNumber.ToString()); 
					scanNumber++;
				}
			}
			return (processedScans); 
		}
	}
}