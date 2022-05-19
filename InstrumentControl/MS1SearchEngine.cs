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
        public short[] ScoreTable;
        protected readonly MS1DatabaseParser Database;
        private int Threads;

        // Constructor will have the basic pieces put together and ready to receive a scan
        public MS1SearchEngine(MS1DatabaseParser database, PpmTolerance tolerance, int threads = 1)
        {
            Tolerance = tolerance;
            Database = database;
            Threads = threads;
        }






        /// <summary>
        /// Will process a MsDataScan and populate the ScoreTable field based upon whether the scan contains proteins within the database based upon Monoisotopic mass
        /// </summary>
        /// <param name="scan">Scan to be processed</param>
        /// <param name="searchType"></param>
        /// <returns>deconvoluted envelopes from the scan passed in</returns>
        public List<IsotopicEnvelope> FindPeakWithinDatabase(MsDataScan scan, string searchType = "boolean")
        {
            // deconvolue scan (eventually switch to unidec || austinDeconV?)
            int minAssumedChargeState = 2;
            int maxAssumedChargeState = 60;
            int deconvolutionTolerancePpm = 6;
            int intensityRatio = 3;
            List<IsotopicEnvelope> envelopes = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).ToList();
            ScoreTable = new short[envelopes.Count];
            List<IsotopicEnvelope> envelopesWithMassInDatabase = new();
            int[] threads = Enumerable.Range(0, Threads).ToArray();

            switch (searchType)
            {
                // score represents the number of masses in the database within tolerance
                case "occurrences":
                    envelopesWithMassInDatabase = envelopes.FindAll(FindInDatabaseWithinMassTolerance);
                    foreach (var match in envelopesWithMassInDatabase)
                    {
                        ScoreTable[envelopes.IndexOf(match)] = (short)Database.ProteinList.Count(p => Tolerance.Within(match.MonoisotopicMass, p.MonoisotopicMass));
                    }
                    break;

                // score of 1 means the mass was found within the database within tolerance
                case "boolean":
                    Parallel.ForEach(threads, index =>
                    {
                        for (; index < envelopes.Count; index += Threads)
                        {
                            if (FindInDatabaseWithinMassTolerance(envelopes[index]))
                            {
                                lock (ScoreTable)
                                {
                                    ScoreTable[index] = 1;
                                }
                            }
                        }
                    });
                    break;
            }
            return envelopes;
        }

        #region Private Helpers

        /// <summary>
        /// Explicit predicate delegate for finding if the envelope exists within the database by monoisotopic mass 
        /// </summary>
        private bool FindInDatabaseWithinMassTolerance(IsotopicEnvelope envelope)
        {
            if (Database.ProteinIndex.Any(p => Tolerance.Within(envelope.MonoisotopicMass, p)))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        private bool FindInDatabaseWithinMassTolerance2(IsotopicEnvelope envelope)
        {
            int index = Array.BinarySearch(Database.ProteinIndex, envelope.MonoisotopicMass);
            if (index >= 0)
            {
                return true;
            }
            else
            {
                index = ~index; // Will be the index of the closest element greater than the search term
                if (index == Database.ProteinIndex.Length && Tolerance.Within(envelope.MonoisotopicMass, Database.ProteinIndex[index - 1]))
                    return true;
                if (index > 0 && index < Database.ProteinIndex.Length)
                {
                    if (Tolerance.Within(envelope.MonoisotopicMass, Database.ProteinIndex[index - 1]))
                        return true;
                    if (Tolerance.Within(envelope.MonoisotopicMass, Database.ProteinIndex[index]))
                        return true;
                }
                return false;
            }
        }

        #endregion
    }
}
