using Chemistry;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;
using Readers; 
using IsdDataProcessing;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestIsdDataProcessing
    {
        [Test]
        public void TestImport()
        {
            string path = @"D:\DeconvolutionPaper\SixProtStandardMixMeth1.raw";
            var scansFull = new Readers.ThermoRawFileReader(path).GetAllScansList();
            
            var itScans = scansFull.GetIonTrapScans();
            // skip the first isfScan. Then interleave. 
            var isfScans = scansFull.GetMs1SidScans().Skip(1);
            var interleaved = itScans.InterleaveScans(isfScans)

        }
    }
}