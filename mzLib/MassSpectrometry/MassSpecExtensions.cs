using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public static class MassSpecExtensions
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="listPeaks"></param>
        public static void LogMzTransformPeaksList(this List<MzPeak> listPeaks)
        {
            for(int i =0; i < listPeaks.Count; i++)
            {
                double mzVal = listPeaks[i].Mz; 
                listPeaks[i].SetNewMz(LogMZTransformedMZ(mzVal)); 
            }
        }
        public static double LogMZTransformedMZ(double mzVal, double adductMass = 1.007)
        {
            return Math.Log(mzVal - 1.007);
        }
    }
}
