using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using NUnit; 
using NUnit.Framework;
using Readers;

namespace Test.TestReaders
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestBruker
    {
        [Test]
        public void TestConstructors()
        {

        }

        [Test]
        public void TestFileDoesntExist()
        {
            string fakePath = "fakePath.d";
            var reader = MsDataFileReader.GetDataFile(fakePath);
            Assert.Throws<FileNotFoundException>(() =>
                reader.InitiateDynamicConnection());
        }

        [Test]
        [TestCase(@"D:\BurkerFileSupport\WT_ESO_+_B_B5_01_4035.d")]
        public void TestLoadAllStaticData(string path)
        {
            MsDataFile brukerData = MsDataFileReader.GetDataFile(path).LoadAllStaticData();
            Assert.That(brukerData.NumSpectra, Is.EqualTo(29767));
            Assert.That(brukerData.Scans[1].Polarity == Polarity.Positive);
            Assert.That(brukerData.Scans[1].DissociationType == DissociationType.CID);
            Assert.That(brukerData.Scans[1].TotalIonCurrent == 4192d);
            Assert.That(brukerData.Scans[1].NativeId == "scan=2");
            Assert.That(brukerData.Scans[1].SelectedIonMZ, Is.EqualTo(426.009674).Within(0.001));
            Assert.That(brukerData.Scans[1].MsnOrder == 2);
            // need to assert that the cumulative total of the y values equals the Total ion current. That would 
            // be a good check for the y values being read correctly. 

        }

        [Test]
        public void TestGetSourceFile()
        {

        }

        [Test]
        public void TestDynamicConnection()
        {

        }

        [Test]
        public void TestPeakFiltering()
        {

        }

        [Test]
        public void TestDynamicRaw()
        {

        }

        [Test]
        public void TestEthcdReading()
        {
            // not sure if I need to do this one. 
        }
    }
}
