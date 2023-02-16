using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media.Animation;
using Accessibility;
using Easy.Common.Extensions;
using Easy.Common.Interfaces;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MathNet.Numerics.Differentiation;
using MathNet.Numerics.Statistics;
using MzLibUtil.MrsNoiseEstimation;
using NUnit.Framework;
using OxyPlot;
using Plotly.NET;
using SpectralAveraging;

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

        [Test]
        public void TestAverageSignalValue()
        {
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        @"D:\Top-Down\HGH\HGH_10uM_2ulmin_StepThroughMeth_20221108164950.RAW")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1
                                && i.MzAnalyzer == MZAnalyzerType.Orbitrap
                                && i.OneBasedScanNumber < 200)
                    .ToList();
            foreach (var scan in scansList)
            {
                MRSNoiseEstimator.MRSNoiseEstimation(scan.MassSpectrum.YArray, 0.01, out double noiseEstimate); 
                double stddevSignal = scan.MassSpectrum.YArray.StandardDeviation();
                double meanOfSignal = scan.MassSpectrum.YArray.Where(i => i > 1.97 * stddevSignal).Average();
                double snr = meanOfSignal / stddevSignal;
            }

            // "Traditionally, SNR is defined to be the ratio of the average signal value to the standard deviation of the signal." 
        }

        [Test]
        public void TestAveragedSNRImprovement()
        {
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        "C:\\Users\\Austin\\source\\repos\\mzLib\\mzLib\\Test\\DataFiles\\small.RAW")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1
                                && i.MzAnalyzer != MZAnalyzerType.Orbitrap)
                    .ToList();
            SpectralAveraging.SpectralAveragingParameters averagingParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01, 
                MaxSigmaValue = 3.0, 

                MinSigmaValue = 0.5, 
                MaxThreadsToUsePerFile = 12, 
                NormalizationType = NormalizationType.AbsoluteToTic, 
                NumberOfScansToAverage = scansList.Count,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping, 
                SpectralWeightingType = SpectraWeightingType.NoiseLevel
            };
            SpectralAveraging.SpectralAveragingParameters noRejectionParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01,
                MaxSigmaValue = 1.5,
                MinSigmaValue = 1.5,
                MaxThreadsToUsePerFile = 12,
                NormalizationType = NormalizationType.AbsoluteToTic,
                SpectralWeightingType = SpectraWeightingType.NoiseLevel,
                NumberOfScansToAverage = scansList.Count,
                OutlierRejectionType = OutlierRejectionType.NoRejection
            };


            var medianRef = scansList[5].MassSpectrum.YArray.Median(); 
            double stdevReference = scansList[5].MassSpectrum.YArray.Where(i => i < medianRef).StandardDeviation();
            double referenceSnr = scansList[5].MassSpectrum.YArray.Where(i => i > stdevReference).Mean() / stdevReference; 
            
            var averagedScan = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(averagingParams);
            var noRejectionAveraged = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(noRejectionParams);

            double averagedScanMedian = averagedScan.YArray.Median(); 
            var averagedStddev = averagedScan.YArray.Where(i => i < averagedScanMedian).StandardDeviation();
            var averagedSnr = averagedScan.YArray.Where(i => i > averagedStddev).Mean() / averagedStddev;

            var rejectionMedian = noRejectionAveraged.YArray.Median(); 
            var noRejectionSD = noRejectionAveraged.YArray.Where(i => i < rejectionMedian).StandardDeviation();
            var noRejectionSnr = noRejectionAveraged.YArray.Where(i => i > noRejectionSD).Mean() / noRejectionSD;

            List<GenericChart.GenericChart> plotList = new();
            GenericChart.GenericChart averagedPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: averagedScan.XArray, y: averagedScan.YArray);
            GenericChart.GenericChart noRejectionPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: noRejectionAveraged.XArray, y: noRejectionAveraged.YArray);
            GenericChart.GenericChart originalPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: scansList[5].MassSpectrum.XArray,
                y: scansList[5].MassSpectrum.YArray);

            plotList.Add(averagedPlot);
            plotList.Add(noRejectionPlot);
            plotList.Add(originalPlot);
            Chart.Combine(plotList).Show();
            Console.WriteLine("Theoretical Improvement = {0}", Math.Sqrt(scansList.Count));

            Console.WriteLine("{0}\t{1}\t{2}", referenceSnr, averagedSnr, noRejectionSnr);
        }
        [Test]
        public void TestAveragedSNRImprovementOrbi()
        {
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        @"D:\Top-Down\KatiesTopDown\02-18-20_jurkat_td_rep2_fract7.raw")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1
                                && i.MzAnalyzer == MZAnalyzerType.Orbitrap
                                && i.RetentionTime > 38 && i.RetentionTime < 38.5) 
                    .ToList();
            SpectralAveraging.SpectralAveragingParameters averagingParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01,
                MaxSigmaValue = 6.0,

                MinSigmaValue = 0.1,
                MaxThreadsToUsePerFile = 1,
                NormalizationType = NormalizationType.AbsoluteToTic,
                NumberOfScansToAverage = scansList.Count,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly
            };
            SpectralAveraging.SpectralAveragingParameters noRejectionParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01,
                MaxSigmaValue = 1.5,
                MinSigmaValue = 1.5,
                MaxThreadsToUsePerFile = 1,
                NormalizationType = NormalizationType.AbsoluteToTic,
                SpectralWeightingType = SpectraWeightingType.WeightEvenly,
                NumberOfScansToAverage = scansList.Count,
                OutlierRejectionType = OutlierRejectionType.NoRejection
            };
            int middleIndex = (scansList.Count - 1) / 2; 

            var medianRef = scansList[middleIndex].MassSpectrum.YArray.Median();
            double stdevReference = scansList[middleIndex].MassSpectrum.YArray.Where(i => i < medianRef).StandardDeviation();
            double referenceSnr = scansList[middleIndex].MassSpectrum.YArray.Where(i => i > stdevReference).Mean() / stdevReference;

            var averagedScan = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(averagingParams);
            var noRejectionAveraged = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(noRejectionParams);

            double averagedScanMedian = averagedScan.YArray.Median();
            double scalingAveraged = medianRef / averagedScanMedian;
            var averagedStddev = averagedScan.YArray
                .Where(i => i < averagedScanMedian)
                .Select(i => i * scalingAveraged)
                .StandardDeviation();
            var averagedSnr = averagedScan.YArray
                .Where(i => i > averagedStddev)
                .Select(i => i * scalingAveraged)
                .Mean() / averagedStddev;

            var rejectionMedian = noRejectionAveraged.YArray.Median();
            var scalingRejection = medianRef / rejectionMedian; 
            var noRejectionSD = noRejectionAveraged.YArray
                .Where(i => i < rejectionMedian)
                .Select(i => i * scalingRejection)
                .StandardDeviation();
            var noRejectionSnr = noRejectionAveraged.YArray
                .Where(i => i > noRejectionSD)
                .Select(i => i * scalingRejection)
                .Mean() / noRejectionSD;

            List<GenericChart.GenericChart> plotList = new();
            GenericChart.GenericChart averagedPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: averagedScan.XArray, y: averagedScan.YArray);
            GenericChart.GenericChart noRejectionPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: noRejectionAveraged.XArray, y: noRejectionAveraged.YArray);
            GenericChart.GenericChart originalPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: scansList[middleIndex].MassSpectrum.XArray,
                y: scansList[5].MassSpectrum.YArray);

            plotList.Add(averagedPlot);
            plotList.Add(noRejectionPlot);
            plotList.Add(originalPlot);
            Chart.Combine(plotList).Show();
            Console.WriteLine("Theoretical Improvement = {0}", Math.Sqrt(scansList.Count));

            Console.WriteLine("{0}\t{1}\t{2}", referenceSnr, averagedSnr, noRejectionSnr);

            Plotly.NET.CSharp.Chart.Histogram<double, double, string>(averagedScan.YArray.Select(i => i * averagedScanMedian / medianRef).ToArray()).Show();
            Plotly.NET.CSharp.Chart.Histogram<double, double, string>(noRejectionAveraged.YArray.Select(i => i * rejectionMedian / medianRef).ToArray()).Show();
            Plotly.NET.CSharp.Chart.Histogram<double, double, string>(scansList[middleIndex].MassSpectrum.YArray)
                .Show(); 

        }
        [Test]
        public void TestAveragedSNRImprovementTof()
        {
            List<MsDataScan> scansList = IO.MzML.Mzml.LoadAllStaticData(
                        @"C:\Users\Austin\source\repos\mzLib\mzLib\Test\DataFiles\WT_ESO_+_B_B5_01_4035.mzML")
                    .GetAllScansList()
                    .Where(i => i.OneBasedScanNumber < 100 
                    && i.MsnOrder == 1)
                    .ToList();
            SpectralAveraging.SpectralAveragingParameters averagingParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01,
                MaxSigmaValue = 3.0,
                MinSigmaValue = 0.5,
                MaxThreadsToUsePerFile = 1,
                NormalizationType = NormalizationType.NoNormalization,
                NumberOfScansToAverage = scansList.Count - 1,
                OutlierRejectionType = OutlierRejectionType.AveragedSigmaClipping,
                SpectralWeightingType = SpectraWeightingType.NoiseLevel
            };
            SpectralAveraging.SpectralAveragingParameters noRejectionParams = new SpectralAveragingParameters()
            {
                BinSize = 0.01,
                MaxSigmaValue = 1.5,
                MinSigmaValue = 1.5,
                MaxThreadsToUsePerFile = 1,
                NormalizationType = NormalizationType.NoNormalization,
                SpectralWeightingType = SpectraWeightingType.NoiseLevel,
                NumberOfScansToAverage = scansList.Count - 1,
                OutlierRejectionType = OutlierRejectionType.NoRejection
            };

            var medianRef = scansList[(scansList.Count - 1)/2].MassSpectrum.YArray.Median();
            double stdevReference = scansList[(scansList.Count - 1) / 2].MassSpectrum.YArray.Where(i => i < medianRef).StandardDeviation();
            double referenceSnr = scansList[(scansList.Count - 1) / 2].MassSpectrum.YArray.Where(i => i > stdevReference).Mean() / stdevReference;

            var averagedScan = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(averagingParams);
            var noRejectionAveraged = scansList.Select(i => i.MassSpectrum).ToList().AverageSpectra(noRejectionParams);

            double averagedScanMedian = averagedScan.YArray.Median();
            var averagedStddev = averagedScan.YArray.Where(i => i < averagedScanMedian).StandardDeviation();
            var averagedSnr = averagedScan.YArray.Where(i => i > averagedStddev).Mean() / averagedStddev;

            var rejectionMedian = noRejectionAveraged.YArray.Median();
            var noRejectionSD = noRejectionAveraged.YArray.Where(i => i < rejectionMedian).StandardDeviation();
            var noRejectionSnr = noRejectionAveraged.YArray.Where(i => i > noRejectionSD).Mean() / noRejectionSD;

            List<GenericChart.GenericChart> plotList = new();
            GenericChart.GenericChart averagedPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: averagedScan.XArray, y: averagedScan.YArray);
            GenericChart.GenericChart noRejectionPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: noRejectionAveraged.XArray, y: noRejectionAveraged.YArray);
            GenericChart.GenericChart originalPlot = Plotly.NET.CSharp.Chart.Line<double, double, string>(x: scansList[(scansList.Count - 1) / 2].MassSpectrum.XArray,
                y: scansList[(scansList.Count - 1) / 2].MassSpectrum.YArray);

            plotList.Add(averagedPlot);
            plotList.Add(noRejectionPlot);
            plotList.Add(originalPlot);
            Chart.Combine(plotList).Show();
            Console.WriteLine("Theoretical Improvement = {0}", Math.Sqrt(scansList.Count));

            Console.WriteLine("{0}\t{1}\t{2}", referenceSnr, averagedSnr, noRejectionSnr);
        }
        /*
         * Spectra are almost indistinguishable in the 
         *
         */
        //[Test]
        //public void PeaksInSpectraArePoisson()
        //{
        //    List<MsDataScan> scansList =
        //        ThermoRawFileReader.LoadAllStaticData(
        //                @"D:\Top-Down\HGH\HGH_10uM_2ulmin_StepThroughMeth_20221108164950.RAW")
        //            .GetAllScansList()
        //            .Where(i => i.MsnOrder == 1
        //                        && i.MzAnalyzer == MZAnalyzerType.Orbitrap
        //                        && i.OneBasedScanNumber < 201)
        //            .ToList();
        //    double mostAbundantMass = 1475.94;
        //    //30,35, 40,45, 50, 60, 70, 80, 90, 100, 200
        //    int[] scanNumberArray = new[] { 2,3,4,5,7,9, 10,11,12,13,14, 15, 20, 25 };

        //    for (int k = 0; k < scanNumberArray.Length; k++)
        //    {
        //        var peaks = scansList.Where(i => i.OneBasedScanNumber - 1 < scanNumberArray[k] + 10)
        //            .Select(i => i.MassSpectrum.YArray[i.MassSpectrum.GetClosestPeakIndex(1475.9439)]).ToList();
        //        var tics = scansList
        //            .Where(i => i.OneBasedScanNumber - 1 < scanNumberArray[k] + 10)
        //            .Select(i => i.TotalIonCurrent).ToList();
        //        var normalizedPeaks = peaks.Select((i, index) => i / tics[index]).ToList();
        //        var peaksNormMean = normalizedPeaks.Mean();
        //        var peaksNormStddev = normalizedPeaks.StandardDeviation();
        //        if (k == scanNumberArray.Length - 1)
        //        {
        //            Console.WriteLine(string.Join("\n", normalizedPeaks));
        //        }
        //        Console.WriteLine(scanNumberArray[k] + " \t" + peaksNormMean + "\t" + peaksNormStddev);
        //    }
        //    // there probably is Poisson noise, but it doesn't matter. 



        //}
        [Test]
        public void CheckSnrUsingIonTrapData()
        {
            double stepSize = 0.001;
            List<MsDataScan> scansList =
                ThermoRawFileReader.LoadAllStaticData(
                        @"C:\\Users\\Austin\\source\\repos\\mzLib\\mzLib\\Test\\DataFiles\\small.RAW")
                    .GetAllScansList()
                    .Where(i => i.MsnOrder == 1
                                && i.MzAnalyzer != MZAnalyzerType.Orbitrap)
                    .ToList();
            Plotly.NET.CSharp.Chart.Histogram<double,double,string>(scansList[0].MassSpectrum.YArray).Show();
            Plotly.NET.CSharp.Chart.Line<double, double, string>(scansList[0].MassSpectrum.XArray, scansList[0].MassSpectrum.YArray).Show();

            double stdev = scansList[0].MassSpectrum.YArray.StandardDeviation(); 
            double mean = scansList[0].MassSpectrum.YArray.Where(i => i > 3.0 * stdev).Average();
            Console.WriteLine(mean / stdev); 



        }

        [Test]
        [TestCase(1)]
        [TestCase(5)]
        [TestCase(10)]
        [TestCase(15)]
        [TestCase(20)]
        public void JurkatNoiseDensity(double ppmTolerance)
        {
            List<MsDataScan> scansList =
            ThermoRawFileReader.LoadAllStaticData(
                @"D:\Top-Down\KatiesTopDown\02-18-20_jurkat_td_rep2_fract7.raw")
            .GetAllScansList()
            .Where(i => i.MsnOrder == 1
                        && i.MzAnalyzer == MZAnalyzerType.Orbitrap 
                        && i.OneBasedScanNumber > 1000 && i.OneBasedScanNumber < 2000)
            .ToList();


            //Console.WriteLine("SharedPeaks\tSpectrum1 Total Peaks\tSpectrum2 Total Peaks\tNon-shared peak density");
            Console.WriteLine("PpmTolerance\tAverageNoiseDensity");
            double noisePeakDensity = 0.0;
            for (int i = 0; i < scansList.Count - 1; i++)
            {
                var indexList = NoiseDensityEstimation.Spectrum1IndexMatches(scansList[i].MassSpectrum,
                    scansList[i + 1].MassSpectrum, ppmTolerance, usePpm: true);
                var countOfSharedPeaks = indexList.Count();
                var totalSpectrum1MzVals = scansList[i].MassSpectrum.XArray.Length;
                var totalSpectrum2MzVals = scansList[i + 1].MassSpectrum.YArray.Length;

                int steps = CalculateTotalTheoreticalMzMeasurements(scansList[i].MassSpectrum.Range.Minimum,
                    scansList[i].MassSpectrum.Range.Maximum, 240000, 400);

                noisePeakDensity += (totalSpectrum1MzVals - countOfSharedPeaks) / (double)steps;
                //Console.WriteLine("{0}\t{1}\t{2}\t{3}",
                //   countOfSharedPeaks, totalSpectrum1MzVals, totalSpectrum2MzVals, (totalSpectrum1MzVals - countOfSharedPeaks) / (double)steps);
            }

            noisePeakDensity /= scansList.Count;
            Console.WriteLine("{0}\t{1}", ppmTolerance, noisePeakDensity);
        }
    }
}
