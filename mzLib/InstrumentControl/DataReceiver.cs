using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Thermo.TNG.Factory;
using Thermo.Interfaces.FusionAccess_V1;
using Thermo.Interfaces.FusionAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.Control.Scans;
using Thermo.Interfaces.SpectrumFormat_V1;
using Thermo.Interfaces.FusionAccess_V1.Control.Scans;
using System.Threading;
using MassSpectrometry; 

namespace InstrumentControl
{
	public class DataReceiver
	{
		Queue<MsDataScan> ScanProcessingQueue = new();
		public void MSScanContainer_MsScanArrived(object sender, MsScanEventArgs e)
		{
			// put scan in queue
			// Keeping IMScan open will hog the memory resources, so you need to get the 
			
			using (IMsScan scan = e.GetScan())
			{
				// need to quickly convert to something else
				// ScanProcessingQueue.Enqueue(scan);
				ScanProcessingQueue.Enqueue(new MsDataScan(scan)); 
			}
			// start new thread for scan in queue (async?)

			// begin processing scan
		}
	}
	public static class IMsScanExtensions
	{
		public static T GetValueFromHeaderDict<T>(this IMsScan imsScan, string headerValue)
		{
			bool success = imsScan.Header.TryGetValue(headerValue, out string value);
			if (success)
			{
				return (T)Convert.ChangeType(value, typeof(T));
			}
			else
			{
				return default(T);
			}

		}
		public static double[,] GetMassSpectrum(this IMsScan scan)
		{
			double[] xarray; // mzs
			double[] yarray; // intensities

			// add error handling for nulls and unequal x and y array length
			// gets the actual data from the IMsScan object
			xarray = scan.Centroids.Select(i => i.Mz).ToArray();
			yarray = scan.Centroids.Select(i => i.Intensity).ToArray();

			// the output array will always have 2 columns. 
			double[,] outputArray = new double[xarray.Length, 2];
			for (int i = 0; i < xarray.Length; i++)
			{
				outputArray[i, 0] = xarray[i];
				outputArray[i, 1] = yarray[i];
			}
			return outputArray;
		}
		public static T GetValueFromTrailerDict<T>(this IMsScan imsScan, string trailerValue)
		{
			bool success = imsScan.Trailer.TryGetValue(trailerValue, out string value);
			if (success)
			{
				return (T)Convert.ChangeType(value, typeof(T));
			}
			else
			{
				return default(T);
			}
		}
	}
}
