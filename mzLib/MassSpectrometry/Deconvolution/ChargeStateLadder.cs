using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public record struct ChargeStateLadder(double Mass, double[] MzVals)
    {
        public double Mass = Mass;
        // guarantees the order is from low m/z to high m/z. 
        public double[] MzVals = MzVals.OrderBy(i => i).ToArray();
    }
}
