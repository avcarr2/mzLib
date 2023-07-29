using Easy.Common.Interfaces;
using MassSpectrometry;

namespace IsdDataProcessing
{
    public static class DataLoading
    {
        public static IEnumerable<MsDataScan> GetIonTrapScans(this List<MsDataScan> scans)
        {
            return scans.Where(i => i.ScanFilter.Contains("ITMS")); 
        }

        public static IEnumerable<MsDataScan> GetMs1SidScans(this List<MsDataScan> scans)
        {
            return scans.Where(i => i.ScanFilter.Contains("sid=100")); 
        }

        public static IEnumerable<MsDataScan> InterleaveScans(this IEnumerable<MsDataScan> ms1s, 
            IEnumerable<MsDataScan> ms2s)
        {
            return ms1s.Zip(ms2s, (f, s) => new[] { f, s }).SelectMany(f => f); 
        }
    }
}