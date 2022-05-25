using IO.MzML;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InstrumentControl
{
    public class TEMPWriteCombinedScansAsmzML
    {


        public static void SaveMergedScanAsMzml(List<MsDataScan> combinedScans, string outputPath)
        {
            SourceFile temp = new SourceFile(null, null, null, null, null);
            MsDataFile combinedScansFile = new MsDataFile(combinedScans.ToArray(), temp);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(combinedScansFile, outputPath, false);

        }


    }
}
