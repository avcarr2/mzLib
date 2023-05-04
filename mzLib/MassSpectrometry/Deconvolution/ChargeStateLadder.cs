using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    internal record struct ChargeStateLadder(double Mass, double[] MzVals)
    {
        public double Mass = Mass;
        public double[] MzVals = MzVals;
    }
}
