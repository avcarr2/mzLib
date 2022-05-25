using System;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using System.Linq;


namespace IMSScanClassExtensions
{
	public static class IMsScanExtensions
	{
		/* imsScan.Header is a <string, string> dictionary. However, some values are
		 * convertible to numerics. Furthermore, there is no guarantee that the keys to the headers values 
		 * will even exist. Therefore, I created an extension to IMsScan based on 
		 * using the IMsScan.Header.TryGetValue. It can be used to implement (in the future) 
		 * error handling in case the header values doesn't exist. 
		 * 
		 * Also, the (T)Convert.ChangeType() portion of the function should hypothetically 
		 * convert the string based on <T>. 
		*/
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

		/*
		 * Function to retrieve the data in the IMsScan object and return a 
		 * double[,] for easy construction into a MzSpectrum object. 
		 */
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

		/*
		 * Does the same thing as GetValueFromHeaderDict, except on the scan trailer 
		 * instead of the header. 
		 */
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
