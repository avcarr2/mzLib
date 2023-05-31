﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class ChargeStateDeconvolutionParams : DeconvolutionParameters
    {
        public int MinCharge { get; }
        public int MaxCharge { get; }
        public double PeakMatchPpmTolerance { get; }
        public double EnvelopeThreshold { get; }
        public double MinimumMassDa { get; set; }
        public int MaxThreads { get; }

        public ChargeStateDeconvolutionParams(int minCharge, int maxCharge, double peakMatchTolerancePpm,
            int maxThreads, double minimumMass = 9500, double envelopeThreshold = 0.6) : base()
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            PeakMatchPpmTolerance = peakMatchTolerancePpm;
            MinimumMassDa = minimumMass;
            MaxThreads = maxThreads;
        }
    }
}