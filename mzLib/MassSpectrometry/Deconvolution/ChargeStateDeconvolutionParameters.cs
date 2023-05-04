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
        public double MinimumMassDa { get; set; }

        public ChargeStateDeconvolutionParams(int minCharge, int maxCharge, double peakMatchTolerancePpm,
            double minimumMass = 9500) : base()
        {
            MinCharge = minCharge;
            MaxCharge = maxCharge;
            PeakMatchPpmTolerance = peakMatchTolerancePpm;
            MinimumMassDa = minimumMass;
        }
    }
}
