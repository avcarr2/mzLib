using System;
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
        public double MinimumMassDa { get; }
        public double MaximumMassDa { get; }
        public int MaxThreads { get; }
        public double DeltaMass { get; }
        public double SequentialChargeStateDiff { get; }
        public double EnvelopeScoreThresh { get; }
        public double PercentageMatchedThresh { get; }



        public ChargeStateDeconvolutionParams(int minCharge, int maxCharge, double peakMatchTolerancePpm,
            int maxThreads, double minimumMass = 9500, double maximumMass = 60000, 
            double envelopeThreshold = 0.6, double deltaMass = 1.003, 
            double sequentialChargeStateDiff = 1d, double envelopeScoreThresh = 0.8, 
            double percentageMatchedThresh = 0d) : base()
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            PeakMatchPpmTolerance = peakMatchTolerancePpm;
            MinimumMassDa = minimumMass;
            MaximumMassDa = maximumMass;
            MaxThreads = maxThreads;
            EnvelopeThreshold = envelopeThreshold;
            DeltaMass = deltaMass;
            SequentialChargeStateDiff = sequentialChargeStateDiff;
            EnvelopeScoreThresh = envelopeScoreThresh;
            PercentageMatchedThresh = percentageMatchedThresh;
        }
    }
}
