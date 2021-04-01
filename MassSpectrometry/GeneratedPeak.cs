﻿using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    internal class GeneratedPeak : MzPeak
    {
        private readonly List<double> mzs = new();
        private readonly List<double> intensities = new();

        public GeneratedPeak(double Mz, double Intensity) : base(Mz, Intensity)
        {
            mzs.Add(Mz);
            intensities.Add(Intensity);
        }

        internal void AddMzPeak(double anotherMz, double anotherIntensity)
        {
            mzs.Add(anotherMz);
            intensities.Add(anotherIntensity);
            Intensity = intensities.Sum();
            double weightedSumMz = 0;
            for (int i = 0; i < mzs.Count; i++)
                weightedSumMz += mzs[i] * intensities[i];
            Mz = weightedSumMz / Intensity;
        }
    }
}