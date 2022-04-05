using System;
using MassSpectrometry;
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
        public void FindPeakWithinDatabase(MsDataScan scan, out List<IsotopicEnvelope> envelopes, out int[] scoreTable, string searchType = "occurrences")
        {
            // deconvolue scan (eventually switch to unidec?)
            int minAssumedChargeState = 2;
            int maxAssumedChargeState = 60;
            int deconvolutionTolerancePpm = 6;
            int intensityRatio = 3;
            envelopes = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).ToList();
            scoreTable = new int[envelopes.Count];
            var envelopesWithMassInDatabase = envelopes.FindAll(FindWithinMass);

            switch (searchType)
            {
                // score represents the number of masses in the database within tolerance
                case "occurrences": 
                    foreach (var match in envelopesWithMassInDatabase)
                    {
                        scoreTable[envelopes.IndexOf(match)] = Database.ProteinList.Count(p => Tolerance.Within(match.MonoisotopicMass, p.MonoisotopicMass));
                    }
                    break;
                // score of 1 means the mass was found within the database within tolerance
                case "boolean":
                    foreach (var match in envelopesWithMassInDatabase)
                    {
                        scoreTable[envelopes.IndexOf(match)] = 1;
                    }
                    break;
            }
        }

        // Explicit predicate delegate for finding if the envelope exists within the database by monoisotopic mass 
        private bool FindWithinMass(IsotopicEnvelope envelope)
        {
            if (Database.ProteinList.Any(p => p.MonoisotopicMass >= Tolerance.GetMinimumValue(envelope.MonoisotopicMass) && p.MonoisotopicMass <= Tolerance.GetMaximumValue(envelope.MonoisotopicMass)))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

    }
}
