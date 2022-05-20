using System;
using NUnit.Framework;
using MassSpectrometry.MzSpectra;
using System.Linq;
using MassSpectrometry; 

namespace Test
{
    public class TestSpectraStandardization
    {
        [Test]
        [TestCase(500.0, 1750.0, 0.001, 1251000)]
        public void TestCreateStandardMZAxis(double lowVal, double highVal, double massAccuraccy, 
            int expectedResult)
        {
            (double,double) range = new (lowVal, highVal);
            double[] resultStdArray = SpectraStandardization.CreateStandardMZAxis(range, massAccuraccy);
            // expected result with this range is 1251000
            Assert.AreEqual(expectedResult, resultStdArray.Length);
            Assert.AreEqual(massAccuraccy, resultStdArray[100] - resultStdArray[99], 0.0005); 
        }
        [Test]
        public void TestLinearInterpolate()
        {
            double x1, y1, x2, y2, x;
            x1 = 500.0;
            y1 = 1.2E6;
            x2 = 500.002;
            y2 = 1E6;
            x = 500.001;

            double result = SpectraStandardization.LinearInterpolate(x1, y1, x2, y2, x);
            Assert.AreEqual(1.1E6, result, 0.01); 
        }
        [Test]
        public void TestCreateNewAxes()
        {
            var function = Enumerable.Range(0, 1000).Select(i => Math.Pow((double)i /100, 0.5)).ToArray();
            // create new array to hold upsample data. 
            double[] resampledData = new double[1200];

            SpectraStandardization.ResampleDataAndInterpolate(function, ref resampledData);
            // Need to assert some things. 
        }
        [Test]
        public void TestCreateStandardScan()
        {
            var stdSpectra = new SpectraStandardization(new(500.0, 1500.0), 0.001);
            double[] testIntVals = Enumerable.Range(0, 1000).Select(i => Math.Pow((double)i / 100, 0.5)).ToArray();
            double[] mzVals = Enumerable.Range(0, 1000).Select(i => (double)i / 100D).ToArray();
            MzSpectrum testSpectrum = new(mzVals, testIntVals, true);
            double[] newXAxis = Enumerable.Range(0, 1200).Select(i => (double)i).ToArray();
            var outputScan = SpectraStandardization.CreateStandardScan(testSpectrum.YArray, newXAxis); 
            Assert.AreEqual(newXAxis.Length, outputScan.YArray.Length); 
        }
    }
}
