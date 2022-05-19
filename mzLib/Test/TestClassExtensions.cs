﻿using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;
using System.IO;
using System.Linq;
using System;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestClassExtensions
    {
        [Test]
        public static void TestBoxCarSmooth()
        {
            double[] inputData = new double[] { 0.19, 0.69, 0.03, 0.85, 0.84, 0.46, 0.09, 0.05, 0.11, 0.5, 0.6, 0.78, 0.48, 0.66, 0.61, 0.78, 0.82, 0.18, 0.77, 0.14, 0.97, 0.48, 0.54, 0.98, 0.01, 0.38, 0.26, 0.4, 0.31, 0.41, 0.03, 0.2, 0.98, 0.36, 0.24, 0.51, 0.14, 0.96, 0.32, 0.9, 0.36, 0.57, 0.97, 0.07, 0.12, 0.73, 0.92, 0.51, 0.04, 0.2, 0.39, 0.32, 0.33, 0.62, 0.32, 0.68, 0.91, 0.3, 0.68, 0.22, 0.89, 0.27, 0.68, 0.08, 0.61, 0.25, 0.82, 0.73, 0.49, 0.76, 0.01, 0.15, 0.13, 0.96, 0.57, 0.58, 0.96, 0.93, 0.5, 0.45, 0.89, 0.44, 0.59, 0.68, 0.71, 0.85, 0.16, 0.18, 0.68, 0.37, 0.22, 0.81, 0.53, 0.26, 0.94, 0.52, 0.66, 0.55, 0.51, 0.14 };
            double[] mySmoothedArray = ClassExtensions.BoxCarSmooth(inputData, 3);
            string[] expectedOutput = new string[] { "0.3", "0.52", "0.57", "0.72", "0.46", "0.2", "0.08", "0.22", "0.4", "0.63", "0.62", "0.64", "0.58", "0.68", "0.74", "0.59", "0.59", "0.36", "0.63", "0.53", "0.66", "0.67", "0.51", "0.46", "0.22", "0.35", "0.32", "0.37", "0.25", "0.21", "0.4", "0.51", "0.53", "0.37", "0.3", "0.54", "0.47", "0.73", "0.53", "0.61", "0.63", "0.54", "0.39", "0.31", "0.59", "0.72", "0.49", "0.25", "0.21", "0.3", "0.35", "0.42", "0.42", "0.54", "0.64", "0.63", "0.63", "0.4", "0.6", "0.46", "0.61", "0.34", "0.46", "0.31", "0.56", "0.6", "0.68", "0.66", "0.42", "0.31", "0.1", "0.41", "0.55", "0.7", "0.7", "0.82", "0.8", "0.63", "0.61", "0.59", "0.64", "0.57", "0.66", "0.75", "0.57", "0.4", "0.34", "0.41", "0.42", "0.47", "0.52", "0.53", "0.58", "0.57", "0.71", "0.58", "0.57", "0.4" };
            string[] actualOutput = mySmoothedArray.Select(v=>Math.Round(v,2).ToString()).ToArray();

            CollectionAssert.AreEquivalent(expectedOutput, actualOutput);
        }
        [Test]
        public void TestFindPeaks()
        {
            double[] vals = Enumerable.Range(0, 1000)
                .Select(i => (double)i/100)
                .ToArray();
            double[] singlePeakFunction = CreateFunction(vals, SinglePeakGeneratorFunction); 
            int[] resultSinglePeak = ClassExtensions.FindPeaks(singlePeakFunction);
            int expectedResults = 500; // ymax is at x = 5 (so index 500) for this function.
            Assert.AreEqual(expectedResults, resultSinglePeak[0]);

            // sine wave with 3 peaks from x = 0 to x = 10.
            double[] sinWavePeaks = CreateFunction(vals, SinWaveGenerator);
            int[] resultingMultiPeak = ClassExtensions.FindPeaks(sinWavePeaks); 
            Assert.AreEqual(4, resultingMultiPeak.Length);

            // cosine wave with four peaks from x = 0 to x = 10. 
            // there is also a peak at x = 0. 
            double[] cosineWavePeaks = CreateFunction(vals, CoSineWaveGenerator);
            int[] resultingCosMultiPeak = ClassExtensions.FindPeaks(cosineWavePeaks);
            Assert.AreEqual(4, resultingCosMultiPeak.Length); 
        }
        public static double[] CreateFunction(double[] values, Func<double, double> function)
        {
            double[] results = new double[values.Length]; 
            for(int i = 0; i < values.Length; i++)
            {
               results[i] = function(values[i]);
            }
            return results; 
        }
        private double SinglePeakGeneratorFunction(double value)
        {
            return 0.3 * Math.Pow(value, 2) + (12.0 * value) - (0.2 * Math.Pow(value, 3));
        }
        public static double SinWaveGenerator(double value)
        {
            return Math.Sin(2 * value); 
        }
        private double CoSineWaveGenerator(double value)
        {
            return Math.Cos(2 * value); 
        }

    }
}