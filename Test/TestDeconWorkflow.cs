﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using UniDecAPI;
using System.IO;
using MassSpectrometry;
using IO.ThermoRawFileReader;
using System.Runtime.InteropServices;

namespace Test
{
	unsafe class TestDeconWorkflow
	{
		private MsDataScan scan;
		private Config config;
		private InputUnsafe inp;
		private float[] xarray;
		private float[] yarray;
		public IntPtr xarrayPtr;
		public IntPtr yarrayPtr;

		[OneTimeSetUp]
		public void Init()
		{
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
			List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			scan = testScan[0];

			// setup the config struct
			UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(scan, out config);

			// setup the input struct
			inp = UniDecAPIMethods.InputMethods.SetupInputs();

			// assign inp the x and y array data
			xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();



			xarrayPtr = Marshal.AllocHGlobal(Marshal.SizeOf(xarray[0]) * xarray.Length);
			Marshal.Copy(xarray, 0, xarrayPtr, xarray.Length); 

			yarrayPtr = Marshal.AllocHGlobal(Marshal.SizeOf(yarray[0]) * yarray.Length);
			Marshal.Copy(yarray, 0, yarrayPtr, yarray.Length); 
			
			inp.dataInt = (float*)yarrayPtr;
			inp.dataMZ = (float*)xarrayPtr;

			inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
			// correctly assigns a value to inp.barr, but the value isn't correct. Need to make sure that it's 
			// char on the C side. 
			UniDecAPIMethods.UtilityMethods.SetLimits(config, inp);  			
		}
		[OneTimeTearDown]
		public void TearDown()
		{
			Marshal.FreeHGlobal(xarrayPtr);
			Marshal.FreeHGlobal(yarrayPtr);
		}

		[Test]
		public void TestUniDecDeconvolutionWorklow()
		{
			int numberOfElementsInBarr = config.lengthmz * config.numz;
			Decon deconResults = new(); 
			byte[] barr = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberOfElementsInBarr);
			/*
			 * inp.barr is a char* pointing to an array of 1s or 0s. These ones and zeroes are stored in an 
			 * 8-bit/1 byte format called ASCII. When converting from C's char to C#'s byte format,
			 * C# is seing the one byte char '49', which is ASCII code for '1' and converting it to 49, instead of 
			 * 1. So we need to convert the ASCII code to the actual number, which means all we need to do to convert
			 * to the correct byte value in C# is subtract the ASCII code for zero, which is 48. The subtraction is implemented in 
			 * the below method. 
			 */
			UniDecAPIMethods.UtilityMethods.ConvertASCIIBytesFromCToByteInCS(ref barr); 

			float threshold = config.psthresh * Math.Abs(config.mzsig) * config.peakshapeinflate;

			int[] starttab = new int[config.lengthmz];
			int[] endtab = new int[config.lengthmz];
			int maxlength = DirectUniDecPort.Convolution.SetStartEnds(config, ref inp, ref starttab, ref endtab, threshold);
			
			int pslen = config.lengthmz * maxlength;
			float[] mzdist = new float[pslen];
			float[] rmzdist = new float[pslen];
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			DirectUniDecPort.MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int[] mind = new int[mlength];
			float[] mdist = new float[mlength]; 

			for(int i = 0; i < mlength; i++)
			{
				mind[i] = i - (mlength - 1) / 2;
				if (config.msig != 0)
				{
					mdist[i] = (float)(Math.Exp(-Math.Pow((i - (zlength - 1) / 2.0), 2)) / (2.0 * config.zsig * config.zsig));
				}
				else
				{
					mdist[i] = 1; 
				}
			}
			int[] zind = new int[zlength];
			float[] zdist = new float[zlength]; 

			int numclose = mlength * zlength;
			int[] closemind = new int[numclose];
			int[] closezind = new int[numclose];
			float[] closeval = new float[numclose];
			int[] closeind = new int[numclose * config.lengthmz * config.numz];
			float[] closearray = new float[numclose * config.lengthmz * config.numz]; 

			for(int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength]; 
			}

			DirectUniDecPort.Normalization.SimpNormSum(mlength, mdist);
			DirectUniDecPort.Normalization.SimpNormSum(zlength, zdist);
			DirectUniDecPort.Normalization.SimpNormSum(numclose, closeval);

			DirectUniDecPort.Blur.MakeSparseBlur(inp, config, numclose, barr, closezind,
				closemind, closeind, closeval, closearray);

			int badness = 1; 
			for(int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if(barr[i] == 1)
				{
					badness = 0; 
				}
			}
			if(badness == 1)
			{
				throw new InvalidOperationException("Badness = 1..."); 
			}

			float dmax = DirectUniDecPort.Blur.MathUtilities.Max(inp.dataInt, config.lengthmz);
			float betafactor = 1;
			if (dmax > 1) 
			{ 
				betafactor = dmax; 
			}
			DirectUniDecPort.FitFunctions.KillB(inp, config, barr);

			deconResults.blur = new float[config.lengthmz * config.numz];
			deconResults.newblur = new float[config.lengthmz * config.numz];
			float[] oldblur = new float[config.lengthmz * config.numz];
			deconResults.baseline = new float[config.lengthmz * config.numz];
			deconResults.noise = new float[config.lengthmz * config.numz];

			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / ((float)(config.numz + 2));
				if (config.baselineflag == 1)
				{

					deconResults.baseline[i] = val;
					deconResults.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (barr[DirectUniDecPort.ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						if (config.isotopemode == 0)
						{
							deconResults.blur[DirectUniDecPort.ArrayIndexing.Index2D(config.numz, i, j)] = val;
						}
						else { deconResults.blur[DirectUniDecPort.ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
					}
					else
					{
						deconResults.blur[DirectUniDecPort.ArrayIndexing.Index2D(config.numz, i, j)] = 0;
					}
				}
			}

			// copy decon.blur to oldblur: 
			oldblur = deconResults.blur;
			// copy decon.newblur to decon.blur: 
			deconResults.newblur = deconResults.blur;

			float[] dataInt2 = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.dataInt, config.lengthmz);

			DirectUniDecPort.Convolution.DeconvolveBaseline(config, inp, deconResults);
			float blurmax = 0F;
			deconResults.conv = 0F;
			int off = 0;

			DirectUniDecPort.Blur.PerformIterations(ref deconResults, config, inp, betafactor, maxlength,
				starttab, endtab, mzdist, numclose, closeind, closearray, zlength, mlength,
				closemind, closezind, mdist, dataInt2, zdist, barr, rmzdist, oldblur);


			if(config.peakshapeinflate != 1 && config.mzsig != 0)
			{
				if(config.speedyflag == 0)
				{
					DirectUniDecPort.MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse: 0, starttab, endtab, maxlength);
				}
				else
				{
					DirectUniDecPort.MZPeak.MakePeakShape1D(config, inp, threshold, mzdist, rmzdist, makeReverse: 0); 
				}
			}

			blurmax = DirectUniDecPort.Blur.MathUtilities.Max(deconResults.blur, config.lengthmz * config.numz);
			float cutoff = 0F; 
			if(blurmax != 0)
			{
				cutoff = 0.000001F; 
			}

			DirectUniDecPort.ArrayIndexing.ApplyCutoff1D(ref deconResults.blur, blurmax * cutoff, config.lengthmz * config.numz);
			// return Decon 
			deconResults.fitdat = new float[config.lengthmz];
			// This line of code is currently not working. 
			deconResults.error = DirectUniDecPort.Blur.ErrorFunctions.ErrFunctionSpeedy(config, deconResults, barr, inp.dataInt, maxlength,
				inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist); 
			if(config.intthresh != -1)
			{
				for(int i = 0; i < config.lengthmz-1; i++)
				{
					if(inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
					{
						deconResults.fitdat[i] = 0F;
						deconResults.fitdat[i + 1] = 0F; 
					}
				}
			}

			if(config.isotopemode == 2)
			{
				DirectUniDecPort.Blur.Isotopes.MonoisotopicToAverageMass(config, inp, deconResults, barr); 
			}

			float newblurmax = blurmax; 
			if(config.rawflag == 0 || config.rawflag == 2)
			{
				if(config.mzsig != 0)
				{
					newblurmax = DirectUniDecPort.Convolution.Reconvolve(config.lengthmz, config.numz, maxlength,
						starttab, endtab, mzdist, deconResults.blur, deconResults.newblur, config.speedyflag, barr); 
				}
			}

		}
	}
}
