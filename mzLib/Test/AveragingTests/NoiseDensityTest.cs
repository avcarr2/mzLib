using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MathNet.Numerics.Differentiation;
using NUnit.Framework;

namespace Test.AveragingTests
{
    public static class NoiseDensityEstimation
    {
        public static IEnumerable<int> CompareMzSpectra(MzSpectrum spectrum1, MzSpectrum spectrum2, double mzMatchTolerance)
        {
            // count number of data points
            int numberPointsSpectrum1 = spectrum1.XArray.Length;
            
            // find indices of m/z values that match within tolerance
            return Spectrum1IndexMatches(spectrum1, spectrum1, mzMatchTolerance); 

        }

        public static IEnumerable<int> Spectrum1IndexMatches(MzSpectrum spectrum1, MzSpectrum spectrum2, double tolerance, 
            Func<double, double, bool> criteriaFunc = null, bool usePpm = true)
        {
            // find first match
            // probably want to iterate over the second array. 
            int spectrum1Indexer = 0;
            int spectrum2Indexer = 0;


            List<int> indicesList = new(); 
            while(spectrum1Indexer < spectrum1.XArray.Length && spectrum2Indexer < spectrum2.XArray.Length)
            {
                
                    double difference = spectrum1.XArray[spectrum1Indexer] - spectrum2.XArray[spectrum2Indexer];
                    if (usePpm)
                    {
                        difference = difference / spectrum1.XArray[spectrum1Indexer] * 1E6; 
                    }

                    if (difference < -tolerance)
                    {
                        // this branch is hit when mz from spectrum 1 is greater than 
                        // mz from spectrum 2. So need to increment the spectrum 2 indexer. 
                        if(++spectrum1Indexer == spectrum1.XArray.Length) break;
                    }else if (difference > tolerance)
                    {
                        // this branch is hit when mz from spectrum 2 is greater than the 
                        // mz from spectrum 1, so need to increment the indexer for spectrum 1. 
                        if (++spectrum2Indexer == spectrum2.XArray.Length) break;
                    }
                    else
                    {
                        if (criteriaFunc == null)
                        {
                            if (spectrum1.YArray[spectrum1Indexer] > 0 && spectrum2.YArray[spectrum2Indexer] > 0)
                            {
                                indicesList.Add(spectrum1Indexer);
                            }
                        }
                        else
                        {
                            if (criteriaFunc.Invoke(spectrum1.YArray[spectrum1Indexer], spectrum2.YArray[spectrum2Indexer]))
                            {
                                indicesList.Add(spectrum1Indexer);
                            }
                        }

                        spectrum1Indexer++;
                        spectrum2Indexer++;
                    }
            }
            return indicesList;
        }

    }

    public class NoiseDensityTest
    {
        [Test]
        public void TestNoiseDensityEstimation()
        {
            double stepSize = 0.001; 
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        @"C:\\Users\\Austin\\source\\repos\\mzLib\\mzLib\\Test\\DataFiles\\small.RAW")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1 
                                && i.MzAnalyzer == MZAnalyzerType.Orbitrap)
                    .ToList(); 

            double maxYValue = scansList[0].MassSpectrum.YArray.Max();

            bool IntensityThresholdPeak(double d, double d1)
            {
                double intensityThreshold = 0.001 * maxYValue;
                return (d >= intensityThreshold || d1 >= intensityThreshold 
                    && d1 / d > 0.8 && d / d1 < 1d / 0.8);
            }

            Console.WriteLine("SharedPeaks\tSpectrum1 Total Peaks\tSpectrum2 Total Peaks\tNon-shared peak density"); 
            for (int i = 0; i < scansList.Count - 1; i++)
            {
                var indexList = NoiseDensityEstimation.Spectrum1IndexMatches(scansList[i].MassSpectrum,
                    scansList[i + 1].MassSpectrum, 0.1, IntensityThresholdPeak);
                var countOfSharedPeaks = indexList.Count();
                var totalSpectrum1MzVals = scansList[i].MassSpectrum.XArray.Length;
                var totalSpectrum2MzVals = scansList[i + 1].MassSpectrum.YArray.Length;
                var rangemz1 = scansList[i].MassSpectrum.Range.Width / stepSize; 
                Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                    countOfSharedPeaks, totalSpectrum1MzVals, totalSpectrum2MzVals, (totalSpectrum1MzVals - countOfSharedPeaks) / 400020d);
            }
        }

        [Test]
        [TestCase(1)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(15)]
        public void TestNoiseDensityEstimateHgh(double ppmTolerance)
        {
            
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        @"D:\Top-Down\HGH\HGH_10uM_2ulmin_StepThroughMeth_20221108164950.RAW")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1
                                && i.MzAnalyzer == MZAnalyzerType.Orbitrap
                                && i.OneBasedScanNumber < 200)
                    .ToList();

            double maxYValue = scansList[0].MassSpectrum.YArray.Max();

            bool IntensityThresholdPeak(double d, double d1)
            {
                double intensityThreshold = 0.0001 * maxYValue;
                return (d >= intensityThreshold || d1 >= intensityThreshold);
            }

            //Console.WriteLine("SharedPeaks\tSpectrum1 Total Peaks\tSpectrum2 Total Peaks\tNon-shared peak density");
            Console.WriteLine("PpmTolerance\tAverageNoiseDensity");
            double noisePeakDensity = 0.0; 
            for (int i = 0; i < scansList.Count - 1; i++)
            {
                var indexList = NoiseDensityEstimation.Spectrum1IndexMatches(scansList[i].MassSpectrum,
                    scansList[i + 1].MassSpectrum, ppmTolerance, IntensityThresholdPeak, usePpm:true);
                var countOfSharedPeaks = indexList.Count();
                var totalSpectrum1MzVals = scansList[i].MassSpectrum.XArray.Length;
                var totalSpectrum2MzVals = scansList[i + 1].MassSpectrum.YArray.Length;
                
                int steps = CalculateTotalTheoreticalMzMeasurements(scansList[i].MassSpectrum.Range.Minimum, 
                    scansList[i].MassSpectrum.Range.Maximum, 120000, 400);

                noisePeakDensity += (totalSpectrum1MzVals - countOfSharedPeaks) / (double)steps;
                //Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                //   countOfSharedPeaks, totalSpectrum1MzVals, totalSpectrum2MzVals, (totalSpectrum1MzVals - countOfSharedPeaks) / (double)steps);
            }

            noisePeakDensity /= scansList.Count; 
            Console.WriteLine("{0}\t{1}", ppmTolerance, noisePeakDensity);
            // add plotting stuff here. 
        }

        [Test]
        public void TestCalculateNumberStepsInMzRange()
        {
            double mz1 = 500d;
            double resolvingPower = 100000;
            double mz = 0;
            double mzPrev = 500d;
            int counter = 0; 
            while (mz <= 2000)
            {
                mz = (mzPrev / resolvingPower / 2) + mzPrev;
                resolvingPower = Math.Sqrt(resolvingPower * resolvingPower * mzPrev / mz);
                mzPrev = mz;
                counter++;
            }
            Console.WriteLine(counter);
        }

        public int CalculateTotalTheoreticalMzMeasurements(double lowMz, double highMz, 
            double resolvingPowerInitial, double mzInitial)
        {
            // update resolving power
            double resolvingPower = Math.Sqrt(resolvingPowerInitial * resolvingPowerInitial * mzInitial / lowMz);
            double mz = 0;
            double mzPrev = lowMz;
            int counter = 0;
            while (mz <= highMz)
            {
                mz = (mzPrev / resolvingPower / 2) + mzPrev;
                resolvingPower = Math.Sqrt(resolvingPower * resolvingPower * mzPrev / mz);
                counter++;
            }
            return counter;
        }
    }
}
