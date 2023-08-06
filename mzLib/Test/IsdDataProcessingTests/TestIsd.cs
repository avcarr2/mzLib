using Chemistry;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;
using Readers; 
using IsdDataProcessing;
using MassSpectrometry;
using UsefulProteomicsDatabases.Generated;

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

        [Test]
        [TestCase(@"D:\AA_08-01-23_ISF\FileList.txt")]
        public void ConvertEntireList(string directoryPath)
        {

            using var streamReader = new StreamReader(directoryPath);

            string? name;
            while ((name = streamReader.ReadLine()) is not null)
            {
                name = name.Trim(); 
                string outputPath = Path.ChangeExtension(name, "mzML");

                var file = new Readers.ThermoRawFileReader(name);
                var scansFull = file.GetAllScansList();


                var itScans = scansFull.GetIonTrapScans();
                // skip the first isfScan. Then interleave. 
                var isfScans = scansFull.GetMs1SidScans().Skip(1);
                var interleaved = itScans.InterleaveScans(isfScans);

                var results = interleaved.UpdateMs2MetaData().UpdateScanStringMetaData().ToArray();

                SourceFile genericSourceFile = new SourceFile("no nativeID format", "mzML format",
                    null, null, null);
                GenericMsDataFile msFile = new GenericMsDataFile(results, genericSourceFile);
                msFile.ExportAsMzML(outputPath, false);
            }


        }
        
    }
}