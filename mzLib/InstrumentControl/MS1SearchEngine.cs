using System;
using MassSpectrometry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;

namespace InstrumentControl
{
    public class MS1SearchEngine
    {
        public PpmTolerance Tolerance { get; }
        protected readonly MS1DatabaseParser Database;

        // Constructor will have the basic pieces put together and ready to receive a scan
        public MS1SearchEngine(MS1DatabaseParser database, PpmTolerance tolerance)
        {
            Tolerance = tolerance;
            Database = database;
        }

        // Will process and the scan and returns a scoreTable with a high number representing a likely chance of the fragment existing in the database
        public void PeakScorer(MsDataScan scan, out List<IsotopicEnvelope> envelopes, out int[] scoreTable)
        {
            // deconvolue scan (eventually switch to unidec?)
            int minAssumedChargeState = 2;
            int maxAssumedChargeState = 60;
            int deconvolutionTolerancePpm = 6;
            int intensityRatio = 3;
            envelopes = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).ToList();

            double maxIntensity = 0;
            if (envelopes.Count() > 0)
            {
                maxIntensity = envelopes.Max(p => p.TotalIntensity);
            }
            scoreTable = new int[envelopes.Count()];

            // Go through scans and increment the score if there is a protein found in database that matches its monomass w/n tolerance
            for (int i = 0; i < envelopes.Count(); i++)
            {
                if (envelopes[i].TotalIntensity <= maxIntensity / 1000) // May need to adjust this number over time
                    continue;
                double mass = envelopes[i].MonoisotopicMass;
                double lowestMassToSearch = Tolerance.GetMinimumValue(mass);
                double highestMassToSearch = Tolerance.GetMaximumValue(mass);

                // TOTRY: Can redo this with a for loop from min to max performing a binary search for the mass in question. Can try after functioning and compare times
                // TOTRY: Instead of using an int, do a simple bool binary -> this will reduce the calcuations needed in the subsequent comparison steps for peak picking but reduces the ability to use other factors in the scoring in the future
                if (Database.ProteinList.Any(p => p.MonoisotopicMass >= lowestMassToSearch && p.MonoisotopicMass <= highestMassToSearch))
                {
                    scoreTable[i] = Database.ProteinList.Count(p => Tolerance.Within(mass, p.MonoisotopicMass));
                }
            }
        }
    }
}
