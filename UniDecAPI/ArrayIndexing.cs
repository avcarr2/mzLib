﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;
using MathNet;

namespace UniDecAPI
{
	// this is extremely bad practice, but I'm pretty desparate at this point. 
	// I'll fix it later. 
	// I think the best path forward is to do all deconvolution steps in C#, but call all the functions 
	// not directly portable from the UniDec DLL. So basically keep all code managed if possible. 

	// TODO: Rearrange for readability. Public methods first. 
	// TODO: Attempt to use managed arrays to pass to the C functions. 
	public partial class DirectUniDecPort
	{
		public static class ArrayIndexing
		{
			public static int Index2D(int ncols, int r, int c)
			{
				return r * ncols + c;
			}
			public static int Index3D(int ncols, int nrows, int r, int c, int d)
			{
				return r * ncols * nrows + c * nrows + d;
			}
			public static int Indexmod(int length, int r, int c)
			{
				int a = c - r;
				int b = length;
				int result = a % b;
				return result < 0 ? result + b : result;
			}
		}
		public static unsafe class Convolution
		{
			public static int SetStartEnds(Config config, ref InputUnsafe inp, ref int[] starttab, ref int[] endtab, float threshold)
			{
				// changed InputUnsafe to pass by reference here, but I think I could've also passed it directly and 
				// created the pointer with a fixed statement as well. Probably better to pass by reference though. 

				int maxlength = 1;
				for (int i = 0; i < config.lengthmz; i++)
				{
					float point = inp.dataMZ[i] - threshold;
					int start, end;
					if (point < inp.dataMZ[0] && config.speedyflag == 0)
					{
						//start = (int)((point - inp->dataMZ[0]) / (inp->dataMZ[1] - inp->dataMZ[0]));
						start = 0 - Nearfast(inp.dataMZ, (float)2 * inp.dataMZ[0] - point, config.lengthmz);
					}
					else
					{
						start = Nearfast(inp.dataMZ, point, config.lengthmz);
					}
					starttab[i] = start;

					point = inp.dataMZ[i] + threshold;
					if (point > inp.dataMZ[config.lengthmz - 1] && config.speedyflag == 0)
					{
						//end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
						end = config.lengthmz - 1 + Nearfast(inp.dataMZ, (float)2 * inp.dataMZ[0] - point, config.lengthmz);
					}
					else
					{
						end = Nearfast(inp.dataMZ, point, config.lengthmz);
					}
					endtab[i] = end;
					if (end - start > maxlength) { maxlength = end - start; }
					//printf("%d %d\n", start, end);
				}
				//printf("Max Length: %d\t", maxlength);
				return maxlength;
			}
			public static int Nearfast(float* dataMz, float point, int numdat)
			{
				return _Nearfast(dataMz, point, numdat);
			}
			public static void DeconvolveBaseline(int lengthmz, float* dataInt, float* baseline, float mzsig)
			{
				_DeconvolveBaseline(lengthmz, dataInt, baseline, mzsig);
			}
			public static float Reconvolve(int lengthmz, int numz, int maxlength, int* starttab, int* endtab,
				float* mzdist, float* blur, float* newblur, int speedyflag, char* barr)
			{
				return _Reconvolve(lengthmz, numz, maxlength, starttab, endtab, mzdist, blur, newblur, speedyflag, barr);
			}
			
			// Private methods 
			[DllImport("TestDLL.dll", EntryPoint = "nearfast")]
			private static extern int _Nearfast(float* dataMz, float point, int numdat);

			[DllImport("TestDLL.dll", EntryPoint = "deconvolve_baseline")]
			private static extern void _DeconvolveBaseline(int lengthmz, float* dataInt, float* baseline,
				float mzsig);
			
			[DllImport("TestDLL.dll", EntryPoint = "Reconvolve")]
			private static extern float _Reconvolve(int lengthmz, int numz, int maxlength, int* starttab, 
				int* endtab, float* mzdist, float* blur, float* newblur, int speedyflag, char* barr);
			
		}
		public static unsafe class MZPeak
		{
			public static void MakePeakShape2D(Config config, InputUnsafe inp, float[] mzdist, float[] rmzdist, 
				 int makereverse, int[] starttab, int[] endtab, int maxlength)
			{
				fixed(float* mzdistPtr = &mzdist[0], rmzdistPtr = &rmzdist[0])
				{
					fixed(int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
					{
						float mzsig = Math.Abs(config.mzsig) * config.peakshapeinflate;
						_MakePeakShape2D(config.lengthmz, maxlength, starttabPtr, endtabPtr,
							inp.dataMZ, mzsig, config.psfun, config.speedyflag, mzdistPtr, rmzdistPtr, makereverse); 
					}
				}
			}
			public static void MakePeakShape1D(float* dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig, int psfun,
				float* mzdist, float* rmzdist, int makereverse)
			{
				_MakePeakShape1D(dataMZ, threshold, lengthmz, speedyflag, mzsig, psfun, mzdist, rmzdist, makereverse);
			}

			[DllImport("TestDLL.dll", EntryPoint = "MakePeakShape2D")]
			private static extern void _MakePeakShape2D(int lengthmz, int maxlength,
				int* starttab, int* endtab, float* dataMZ, float mzsig, int psfun, int speedyflag,
				float* mzdist, float* rmzdist, int makereverse);

			[DllImport("TestDLL.dll", EntryPoint = "MakePeakShape1D")]
			private static extern void _MakePeakShape1D(float* dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig,
				int psfun, float* mzdist, float* rmzdist, int makereverse); 
			
		}
		public unsafe static class Normalization
		{
			public static void SimpNormSum(int length, float[] data)
			{
				fixed(float* dataPtr = &data[0])
				{
					simp_norm_sum(length, dataPtr);
				}
			}

			[DllImport("TestDLL.dll", EntryPoint = "simp_norm_sum")]
			private static extern void simp_norm_sum(int length, float* data);

		}
		public static unsafe class Blur
		{
			public static void MakeSparseBlur(int numclose, byte* barr, int* closezind,
				int* closemind, float* mtab, int* nztab, float* dataMZ, int* closeind, float* closeval,
				float* closearray, Config config)
			{
				_MakeSparseBlur(numclose, barr, closezind, closemind, mtab,
					nztab, dataMZ, closeind, closeval, closearray, config);
			}
			public static void MakeSparseBlur(InputUnsafe inp, Config config, int numclose, byte[] barr, 
				int[] closezind, int[] closemind, int[] closeind, 
				float[] closeval, float[] closearray)
			{
				fixed(byte* barrPtr = &barr[0])
				{
					fixed(int* closezindPtr = &closezind[0], closemindPtr = &closemind[0], 
						closeindPtr = &closeind[0])
					{
						fixed(float* closevalPtr = &closeval[0], closearrayPtr = &closearray[0])
						{
							_MakeSparseBlur(numclose, barrPtr, closezindPtr, closemindPtr, inp.mtab, inp.nztab,
								inp.dataMZ, closeindPtr, closevalPtr, closearrayPtr, config); 
						}
					}
				}
			}
			public static void CreateInitialBlur(Decon decon, InputUnsafe inp, Config config)
			{
				for (int i = 0; i < config.lengthmz; i++)
				{
					float val = inp.dataInt[i] / (float)(config.numz + 2);
					if (config.baselineflag == 1)
					{

						decon.baseline[i] = val;
						decon.noise[i] = val;
					}

					for (int j = 0; j < config.numz; j++)
					{
						if (inp.barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
						{
							if (config.isotopemode == 0)
							{
								decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = val;
							}
							else { decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
						}
						else
						{
							decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 0;
						}
					}
				}
			}
			public static void SoftargmaxTransposed(float* blur, int lengthmz, int numz, float beta,
				byte* barr, int maxlength, int isolength, int* isotopepos, float* isotopeval, int speedyflag,
				int* starttab, int* endtab, float* mzdist, float mzsig)
			{
				_SoftArgmaxTransposed(blur, lengthmz, numz, beta, barr, maxlength,
					isolength, isotopepos, isotopeval, speedyflag, starttab, endtab, mzdist, mzsig);
			}
			public static void Softargmax(float* blur, int lengthmz, int numz, float beta)
			{
				_Softargmax(blur, lengthmz, numz, beta);
			}
			public static void PointSmoothing(float* blur, byte* barr, int lengthmz, int numz, int width)
			{
				_PointSmoothing(blur, barr, lengthmz, numz, width);
			}
			public static void PointSmoothingPeakWidth(int lengthmz, int numz, int maxlength, int* starttab, int* endtab, float* mzdist, float* blur, int speedyflag, byte* barr)
			{
				_PointSmoothingPeakWidth(lengthmz, numz, maxlength, starttab, endtab, mzdist, blur, speedyflag, barr);
			}
			public static void BlurItMean(int lengthmz, int numz, int numclose, int* closeind, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog)
			{
				_BlurItMean(lengthmz, numz, numclose, closeind, newblur,
				blur, barr, closearray, zerolog);
			}
			public static void BlurItHybrid1(int lengthmz, int numz, int zlength, int mlength,
				int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog)
			{
				_BlurItHybrid1(lengthmz, numz, zlength, mlength, closeind, closemind, closezind,
					mdist, zdist, newblur, blur, barr, closearray, zerolog);
			}
			public static void BlurIt(int lengthmz, int numz, int numclose, int* closeind,
				float* closearray, float* newblur, float* blur, byte* barr)
			{
				_BlurIt(lengthmz, numz, numclose, closeind, closearray, newblur, blur, barr);
			}
			public static void PerformIterations(Decon decon, Config config, InputUnsafe inp, float betafactor, int maxlength,
				int* starttab, int* endtab, float* mzdist, int numclose, int* closeind, float* closearray, int zlength, int mlength,
				int* closemind, int* closezind, float* mdist, float* dataInt2, float* zdist, byte* barr, float* rmzdist, float* oldblur)
			{
				int off = 0;
				for (int iterations = 0; iterations < Math.Abs(config.numit); iterations++)
				{
					decon.iterations = iterations;
					if (config.beta > 0 && iterations > 0)
					{

						Softargmax(decon.blur, config.lengthmz, config.numz, config.beta / betafactor);
						//printf("Beta %f\n", beta);
					}
					else if (config.beta < 0 && iterations > 0)
					{
						SoftargmaxTransposed(decon.blur, config.lengthmz, config.numz, Math.Abs(config.beta / betafactor), inp.barr, maxlength, config.isolength, inp.isotopeops, inp.isotopeval, config.speedyflag, starttab, endtab, rmzdist, config.mzsig);
					}

					if (config.psig >= 1 && iterations > 0)
					{
						PointSmoothing(decon.blur, inp.barr, config.lengthmz, config.numz, Math.Abs((int)config.psig));
						//printf("Point Smoothed %f\n", config.psig);
					}
					else if (config.psig < 0 && iterations > 0)
					{
						PointSmoothingPeakWidth(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, config.speedyflag, barr);
					}


					//Run Blurs
					if (config.zsig >= 0 && config.msig >= 0)
					{
						BlurItMean(config.lengthmz, config.numz, numclose, closeind, decon.newblur, decon.blur, inp.barr, closearray, config.zerolog);
					}
					else if (config.zsig > 0 && config.msig < 0)
					{
						BlurItHybrid1(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
					}
					else if (config.zsig < 0 && config.msig > 0)
					{
						BlurItHybrid2(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
					}
					else
					{
						BlurIt(config.lengthmz, config.numz, numclose, closeind, closearray, decon.newblur, decon.blur, barr);
					}

					//Run Richardson-Lucy Deconvolution
					DeconvolveIterationSpeedy(config.lengthmz, config.numz, maxlength,
						decon.newblur, decon.blur, barr, config.aggressiveflag, dataInt2,
						config.isolength, inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist, rmzdist, config.speedyflag,
						config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);

					//Determine the metrics for conversion. Only do this every 10% to speed up.
					if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations > 0.9 * config.numit))
					{
						float diff = 0;
						float tot = 0;
						for (int i = 0; i < config.lengthmz * config.numz; i++)
						{
							if (barr[i] == 1)
							{
								diff += (float)Math.Pow(((double)decon.blur[i] - (double)oldblur[i]), 2);
								tot += decon.blur[i];
							}
						}
						if (tot != 0) { decon.conv = (diff / tot); }
						else { decon.conv = 12345678; }

						//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
						if (decon.conv < 0.000001)
						{
							if (off == 1 && config.numit > 0)
							{

								break;
							}
							off = 1;
						}
					}
				}
			}
			[DllImport("TestDLL.dll", EntryPoint = "MakeSparseBlur")]
			private static extern void _MakeSparseBlur(int numclose, byte* barr, int* closezind,
				int* closemind, float* mtab, int* nztab, float* dataMZ, int* closeind, float* closeval,
				float* closearray, Config config);

			[DllImport("TestDLL.dll", EntryPoint = "softargmax")]
			private static extern void _Softargmax(float* blur, int lengthmz, int numz, float beta);

			[DllImport("TestDLL.dll", EntryPoint = "softargmax_transposed")]
			private static extern void _SoftArgmaxTransposed(float* blur, int lengthmz, int numz, float beta,
				byte* barr, int maxlength, int isolength, int* isotopepos, float* isotopeval, int speedyflag,
				int* starttab, int* endtab, float* mzdist, float mzsig);

			[DllImport("TestDLL.dll", EntryPoint = "point_smoothing")]
			private static extern void _PointSmoothing(float* blur, byte* barr, int lengthmz, int numz, int width);

			[DllImport("TestDLL.dll", EntryPoint = "point_smoothing_peak_width")]
			private static extern void _PointSmoothingPeakWidth(int lengthmz, int numz, int maxlength, int* starttab, int* endtab, float* mzdist, float* blur, int speedyflag, byte* barr);

			[DllImport("TestDLL.dll", EntryPoint = "blur_it_mean")]
			private static extern void _BlurItMean(int lengthmz, int numz, int numclose, int* closeind, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog);

			[DllImport("TestDLL.dll", EntryPoint = "blur_it_hybrid1")]
			private static extern void _BlurItHybrid1(int lengthmz, int numz, int zlength, int mlength,
				int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog);
			
			[DllImport("TestDLL.dll", EntryPoint = "blur_it_hybrid2")]
			private static extern void _BlurItHybrid2(int lengthmz, int numz, int zlength, int mlength,
				int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog);
			public static void BlurItHybrid2(int lengthmz, int numz, int zlength, int mlength,
				int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
				float* blur, byte* barr, float* closearray, float zerolog)
			{
				_BlurItHybrid2(lengthmz, numz, zlength, mlength, closeind, closemind, closezind,
					mdist, zdist, newblur, blur, barr, closearray, zerolog);
			}
			[DllImport("TestDLL.dll", EntryPoint = "blur_it")]
			private static extern void _BlurIt(int lengthmz, int numz, int numclose, int* closeind,
				float* closearray, float* newblur, float* blur, byte* barr);
			
			[DllImport("TestDLL.dll", EntryPoint = "deconvolve_iteration_speedy")]
			private static extern void _DeconvolveIterationSpeedy(int lengthmz, int numz, int maxlength,
				float* blur, float* blur2, byte* barr, int aggressiveflag, float* dataInt,
				int isolength, int* isotopepos, float* isotopeval, int* starttab, int* endtab,
				float* mzdist, float* rmzdist, int speedyflag, int baselineflag, float* baseline,
				float* noise, float mzsig, float* dataMZ, float filterwidth, float psig);
			public static void DeconvolveIterationSpeedy(int lengthmz, int numz, int maxlength,
				float* blur, float* blur2, byte* barr, int aggressiveflag, float* dataInt,
				int isolength, int* isotopepos, float* isotopeval, int* starttab, int* endtab,
				float* mzdist, float* rmzdist, int speedyflag, int baselineflag, float* baseline,
				float* noise, float mzsig, float* dataMZ, float filterwidth, float psig)
			{
				_DeconvolveIterationSpeedy(lengthmz, numz, maxlength,
				blur, blur2, barr, aggressiveflag, dataInt, isolength, isotopepos, isotopeval,
				starttab, endtab, mzdist, rmzdist, speedyflag, baselineflag, baseline,
				noise, mzsig, dataMZ, filterwidth, psig);
			}

			public static unsafe class MathUtilities
			{
				[DllImport("TestDLL.dll", EntryPoint = "Max")]
				private static extern float _Max(float* blur, int length);
				public static float Max(float* blur, int length)
				{
					return _Max(blur, length);
				}
				public static void AplyCutoff1D(float* array, float cutoff, int lengthmz)
				{
					for(int i = 0; i < lengthmz; i++)
					{
						if (array[i] < cutoff)
						{
							array[i] = 0; 
						}
					}
				}
			}
			public static unsafe class Isotopes
			{
				[DllImport("TestDLL.dll", EntryPoint = "setup_and_make_isotopes")]
				private static extern void _SetupAndMakeIsotopes(ref Config config, ref InputUnsafe inp);
				public static void SetupAndMakeIsotopes(ref Config config, ref InputUnsafe inp)
				{
					_SetupAndMakeIsotopes(ref config, ref inp);
				}
				public static void MonoisotopicToAverageMass(int lengthmz, int numz, float* blur, char* barr,
					int isolength, int* isotopepos, float* isotopeval)
				{
					float[] newblur = new float[lengthmz * numz];
					for (int i = 0; i < lengthmz; i++)
					{
						for (int j = 0; j < numz; j++)
						{
							if (barr[ArrayIndexing.Index2D(numz, i, j)] == 1)
							{
								float topval = blur[ArrayIndexing.Index2D(numz, i, j)];
								for (int k = 0; k < isolength; k++)
								{
									int pos = isotopepos[ArrayIndexing.Index3D(numz, isolength, i, j, k)];
									float val = isotopeval[ArrayIndexing.Index3D(numz, isolength, i, j, k)];
									newblur[ArrayIndexing.Index2D(numz, pos, j)] += topval * val;

								}
							}
						}
					}
				}
			}
			public static unsafe class ErrorFunctions
			{
				[DllImport("TestDLL.dll", EntryPoint = "errfunspeedy")]
				private static extern float _ErrFunctionSpeedy(Config config, Decon decon, char* barr, float* dataInt, int maxlength,
					int* isotopepos, float* isotopeval, int* starttab, int* endtab, float* mzdist, float* rsquared); 
				public static float ErrFunctionSpeedy(Config config, Decon decon, char* barr, float* dataInt, int maxlength,
					int* isotopepos, float* isotopeval, int* starttab, int* endtab, float* mzdist, float* rsquared)
				{
					return _ErrFunctionSpeedy(config, decon, barr, dataInt, maxlength, isotopepos, isotopeval, 
						starttab, endtab, mzdist, rsquared); 
				}
			}
			public static unsafe class MassIntensityDetermination
			{
				public static void IntegrateMassIntensities(Config config, Decon decon, InputUnsafe inp)
				{
					float massmax = config.masslb;
					float massmin = config.massub; 
					if (config.poolflag == 0)
					{
						if (config.rawflag == 1 || config.rawflag == 3)
						{
							IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.blur, decon.massgrid);
						}
						if (config.rawflag == 0 || config.rawflag == 2)
						{
							IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.newblur, decon.massgrid);
						}
					}
					else if (config.poolflag == 1)
					{
						if (config.rawflag == 1 || config.rawflag == 3)
						{
							InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
								inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
						}
						if (config.rawflag == 0 || config.rawflag == 2)
						{
							InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
								inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
						}
					}
					else if (config.poolflag == 2)
					{
						if (config.rawflag == 1 || config.rawflag == 3)
						{
							SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
								inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
						}
						if (config.rawflag == 0 || config.rawflag == 2)
						{
							SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
								inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
						}
					}
					else
					{
						// throw new error
					}
				}
				[DllImport("TestDLL.dll", EntryPoint = "IntegrateTransform")]
				private static extern void _IntegrateTransform(int lengthmz, int numz, float* mtab, float massmax, float massmin, int maaxle, 
					float* massaxis, float* massaxisval, float* blur, float* massgrid);
				public static void IntegrateTransform(int lengthmz, int numz, float* mtab, float massmax, float massmin, int maaxle,
					float* massaxis, float* massaxisval, float* blur, float* massgrid)
				{
					_IntegrateTransform(lengthmz, numz, mtab, massmax, massmin, maaxle,
					massaxis, massaxisval, blur, massgrid); 
				}
				[DllImport("TestDLL.dll", EntryPoint = "InterpolateTransform")]
				private static extern void _InterpolateTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
					float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur);
				public static void InterpolateTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
					float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur)
				{
					_InterpolateTransform(maaxle, numz, lengthmz, nztab, massaxis, adductmass, 
						dataMZ, massgrid, massaxisval, blur); 
				}
				[DllImport("TestDLL.dll", EntryPoint = "SmartTransform")]
				private static extern void _SmartTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
					float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur);
				public static void SmartTransform(int maaxle, int numz, int lengthmz, int* nztab, float* massaxis,
					float adductmass, float* dataMZ, float* massgrid, float* massaxisval, float* blur)
				{
					_SmartTransform(maaxle, numz, lengthmz, nztab, massaxis, adductmass, 
						dataMZ, massgrid, massaxisval, blur); 
				}
			}
			public static unsafe class Scoring
			{
				[DllImport("TestDLL.dll", EntryPoint = "score")]
				private static extern float _UniScore(Config config, Decon decon, InputUnsafe inp, float scoreThreshold); 
				public static float UniScore(Config config, Decon decon, InputUnsafe inp, float scoreThreshold)
				{
					return _UniScore(config, decon, inp, scoreThreshold); 
				}
			}
		}
	}
}
