using Chemistry;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;
using Readers; 
using IsdDataProcessing;
using MassSpectrometry;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestIsdDataProcessing
    {
        [Test]
        public void TestImportExport()
        {
            string path = @"D:\DeconvolutionPaper\SixProtStandardMixMeth1.raw";
            string path2 = @"D:\DeconvolutionPaper\SixProtStandardMixMeth1.mzML";

            var file = new Readers.ThermoRawFileReader(path);
            var scansFull = file.GetAllScansList(); 


            var itScans = scansFull.GetIonTrapScans();
            // skip the first isfScan. Then interleave. 
            var isfScans = scansFull.GetMs1SidScans().Skip(1);
            var interleaved = itScans.InterleaveScans(isfScans);

            var results = interleaved.UpdateMs2MetaData().UpdateScanStringMetaData().ToArray();

            SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                null, null, null);
            GenericMsDataFile msFile = new GenericMsDataFile(results, genericSourceFile);
            msFile.ExportAsMzML(path2, false);
        }
    }
}