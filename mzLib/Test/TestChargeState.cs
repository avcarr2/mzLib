// Copyright 2012, 2013, 2014 Derek J. Bailey
//
// This file (MassTestFixture.cs) is part of CSMSL.Tests.
//
// CSMSL.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// CSMSL.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with CSMSL.Tests. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using NUnit.Framework;
using System;
using Proteomics;
using MassSpectrometry;
using System.IO;
using IO.ThermoRawFileReader;
using MassSpectrometry.MzSpectra;
using System.Collections.Generic;
using System.Linq;
using IO.MzML;
using IO.MzML; 


namespace Test
{
    public class TestChargeState
    {
        public MzSpectrum Scan { get; set; } 
        public LogTransformedSpectra LogMzScan { get; set; }

        [OneTimeSetUp]
        public void CreateScan()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon_centroided.mzML");
            Scan = Mzml.LoadAllStaticData(path).GetAllScansList()[0].MassSpectrum;
            LogMzScan = new LogTransformedSpectra(Scan, 1.007); 
        }         

        [Test]
        public void TestScanForChargeStates()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
            MzSpectrum scan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList()[0].MassSpectrum;

        }
        [Test]
        public void ExampleChargeStateSelection()
		{   
            /*
             * m/z+1 will always be a smaller value than m/z. 
             * So, by using Math.Truncate(m/z - m/z+1), we get the z value of m/z, 
             * and the next integer value up will be the value of the z + 1 charge state. 
             */
            double peak1 = Math.Log(1318.28) + Math.Log(1.007);
            double peak2 = Math.Log(1153.62) + Math.Log(1.007);

            double difference = peak1 - peak2; 
            double recipDiff = -1 / difference;
            Console.WriteLine(recipDiff.ToString()); 
        }
        [Test]
        public void ExampleChargeStateSelectionFromOriginalChargeState()
        {
            /*
             * m/z+1 will always be a smaller value than m/z. 
             * So, by using Math.Truncate(m/z - m/z+1), we get the z value of m/z, 
             * and the next integer value up will be the value of the z + 1 charge state. 
             */
            double peak1 = Math.Log(1318.28) + Math.Log(1.007);
            double peak2 = Math.Log(1153.62) + Math.Log(1.007);

            double difference = peak1 - peak2;
            double recipDiff = 1 / difference;

            int initialZ = (int)Math.Truncate(recipDiff);
            double initialZPlusOne = initialZ + 1;

            // calculate the range of log(m/z) values to look for peaks. 
            // the next (lower) log(m/z) value will be log(z2/z1) 
            double initialZPlusTwo = initialZPlusOne + 1;
            double nextLogMz = Math.Log(initialZPlusTwo / initialZPlusOne);
            Console.WriteLine(nextLogMz);
            double peak3predicted = peak2 - nextLogMz;
            Console.WriteLine(Math.Exp(peak3predicted - Math.Log(1.007)));

            // the next higher log(m/z) value will be log(z1/z2)
            double initialZMinus1 = initialZ - 1;
            double nextHigherLogMz = Math.Log(initialZ / initialZMinus1);
            double peaksMinus1Predicted = peak1 + nextHigherLogMz;
            Console.WriteLine(Math.Exp(peaksMinus1Predicted - Math.Log(1.007))); 

            // The goal is to iterate through the list of peaks in a spectra, assigning each peak to a particular 
            // charge state. 

        }

        // IsotopeCluster
        [Test]
        public void TestGetAllPairs()
        {
            List<int> pairs = new(); 
            pairs.Add(1);
            pairs.Add(2);   
            pairs.Add(3);

            List<(int,int)> pairsList = IsotopeCluster.GetAllPairs(pairs).ToList();
            Assert.AreEqual(3, pairsList.Count);
            Assert.AreEqual((2,3), pairsList[2]); 
        }
        [Test]
        public void TestIsotopomerCalculation()
        {
            double logMzDiff = 0.0008339;
            double adductMass = 1.007;
            double isotopeMass = 1.003; 
            double result = IsotopeCluster.IsotopomerCalculation(logMzDiff, 
                isotopeMass, adductMass);
            double expected = 1.0030;
            Assert.AreEqual(expected, result, 0.001); 
        }
        [Test]
        public void TestPeakIsIsotopeMassAway()
        {
            double peak1 = Math.Log(401.0070) + Math.Log(1.007);
            double peak2 = Math.Log(401.3413) + Math.Log(1.007);

            bool testBool = IsotopeCluster.PeakIsIsotopeMassAway(peak1, peak2,
                1.00335, 0.001, 1.007);
            Assert.AreEqual(true, testBool); 
        }
        [Test]
        public void TestIsotopeCluster()
        {
            List<MzPeak> peaks = new();
            peaks.Add(new MzPeak(401.0070, 1E6));
            peaks.Add(new MzPeak(401.3413, 1.2E6));
            peaks.Add(new MzPeak(401.6757, 1.0E6));
            var logMzPeaks = peaks.Select(i => new MzPeak(Math.Log(i.Mz) + Math.Log(1.007),i.Intensity))
            .ToList(); 

            IsotopeCluster cluster = new(logMzPeaks);
            Assert.AreEqual(3, cluster.PeaksInClusterList.Count);
            Assert.AreEqual(1.2E6, cluster.MaxPeakInCluster.Intensity);
            Assert.AreEqual(3, cluster.NumberPeaksInCluster);
            Assert.AreEqual(3, cluster.ChargeState); 
        }
        [Test]
        public void TestAddProposedPeaksAndFilterValidity()
        {
            List<MzPeak> peaks = new();
            peaks.Add(new MzPeak(401.0070, 1E6));
            peaks.Add(new MzPeak(401.3413, 1.2E6));
            peaks.Add(new MzPeak(401.6757, 1.0E6));
            // next 2 are peaks that are not supposed to belong: 
            peaks.Add(new MzPeak(400.000, 2E4));
            peaks.Add(new MzPeak(401.30, 5E5)); 
            var logMzPeaks = peaks.Select(i => new MzPeak(Math.Log(i.Mz) + Math.Log(1.007), i.Intensity))
            .ToList();

            IsotopeCluster cluster = new();
            cluster.AddProposedPeaksAndFilterValidity(logMzPeaks); 
            Assert.AreEqual(3, cluster.PeaksInClusterList.Count);
            Assert.AreEqual(1.2E6, cluster.MaxPeakInCluster.Intensity);
            Assert.AreEqual(3, cluster.NumberPeaksInCluster);
            Assert.AreEqual(3, cluster.ChargeState);
        }
        [Test]
        [TestCase(0.01)]
        [TestCase(0.1)]
        [TestCase(0.05)]
        [TestCase(0.5)]
        [TestCase(1.0)]
        [TestCase(1.1)]
        [TestCase(1.2)]
        public void TestCheckValidChargeState(double stdev)
        {
            double[] chargeStateGuesses = new double[] {4.44, 6.776, 6.776,
                7.09, 7.09, 7.09, 6.77, 7.09, 7.02, 6.77, 7.08, 7.08, 4.708,
                9.13, 8.611};
            IsotopeCluster.CheckValidChargeState(chargeStateGuesses, stdev, out int charge);
            Assert.AreEqual(7, charge); 
        }
        [Test]
        public void TestCalculateChargeState()
        {

        }
        // LogTransformedSpectra
        [Test]
        public void Test()
        {

        }
        // ChargeStateEnvelope
        [Test]
        [TestCase(0.001)]
        //[TestCase(0.0025)]
        //[TestCase(0.01)]
        public void TestFindInitialIsotopicCluster(double window)
        {
            ChargeStateEnvelope cse = new();
            // range is too wide and there's still too many false positive 
            // isotopomer peaks.  
            cse.FindInitialIsotopicCluster(LogMzScan, window);
            // currently fails to get the correct charge of the isotopic envelope, 
            // which is 7, because the centroiding has slight shifts between peaks. 
            Console.WriteLine(cse.IsotopeClusters[0].NumberPeaksInCluster); 
        }
        
        [Test]
        public void TestHelpAshley()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "mNP_DIA_DSD_13.mzML");
            var scans = Mzml.LoadAllStaticData(path).GetAllScansList();
            // list of all the ion injection times for each ms1 scan

            var injectionTimes = scans
                .Where(i => i.MsnOrder == 1)
                .Select(i => i.InjectionTime);
            Console.WriteLine(String.Join("\n", injectionTimes));
        }
    }
}