using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class AustinConv
    {
        public int[] Charges { get; set; }
        public double[] LogZ { get; set; }
        // driver code 
        public void RunDeconvolution(int minCharge, int maxCharge)
        {
            CreatesChargesAndLogZ(minCharge, maxCharge); 
            // smooth spectrum 
            // find peaks from smoothed spectrum
            // find all putatitive charge state envelopes
            // Add the real peaks to the charge state envelopes
            List<ChargeStateEnvelope> cse = new List<ChargeStateEnvelope>();
            // create the isotopic envelopes 
        }
        public void FindAllChargeStateEnvelopes(List<MzPeak> peaksList, double zeroThreshold, double peakFindThresh, 
            int maxIterations)
        {
            List<List<(int, double)>> cseList = new();
            int iterations = 0; 
            while(peaksList.Count > 0 && iterations < maxIterations) 
            {
                var chargeMZPairs = FindPutativeChargeStateEnvelope(peaksList, zeroThreshold, peakFindThresh);
                cseList.Add(chargeMZPairs);
                UpdatePeakList(peaksList, chargeMZPairs.Select(i => i.Item2).ToList()); 
                iterations++;
            }

        }
        /// <summary>
        /// Takes list of MzPeak representing the peaks from the smoothed spectrum and 
        /// processes them into list of related charge states. 
        /// </summary>
        /// <param name="peaksList"></param>
        public List<(int, double)> FindPutativeChargeStateEnvelope(List<MzPeak> peaksList, double zeroThreshold, double peakFindThreshold)
        {
            List<(int, double)> chargeMzPairs = new(); 
            MzPeak maxPeak = peaksList.MaxBy(i => i.Intensity);
            // find initial charge state pair
            bool csFound = FindInitialChargeStatePair(maxPeak, peaksList, zeroThreshold, out (int, double)? nextPeak);
            // if found, find all other peaks that correspond to that charge state
            if (csFound)
            {
                FindAllRelatedChargeStates(nextPeak.Value, peaksList, peakFindThreshold, chargeMzPairs); 
            }
            return chargeMzPairs; 
        }
        public bool FindInitialChargeStatePair(MzPeak maxPeak, List<MzPeak> peaksList, 
            double zeroThreshold, out (int,double)? nextPeak)
        {
            // max peak distance is the log(minCharge+1/minCharge)
            // only search peaks lower in m/z than the max peak 
            double maxPeakDistance = LogZ[1] - LogZ[0];

            double maxPeakDistToInclude = maxPeak.Mz - maxPeakDistance;
            // filter peaks 
            var filteredPeaksList = peaksList
                .Where(i => i.Mz >= maxPeakDistance && i.Mz < maxPeak.Mz)
                .ToList();
            List<(int, double)> validCharges = new();
            for (int i = 0; i < filteredPeaksList.Count; i++)
            {
                // search for the correct charge state would probably be improved 
                // by a binary search-like procedure

                for (int j = 0; j < Charges.Length; j++)
                {
                    bool success = IsConsecutiveChargeState(maxPeak.Mz, filteredPeaksList[i].Mz, 
                        Charges[j], zeroThreshold);
                    if (success)
                    {
                        validCharges.Add((Charges[j], filteredPeaksList[i].Mz)); 
                    }
                }
            }    
            if(validCharges.Count > 0)
            {
                if(validCharges.Count == 1)
                {
                    nextPeak = validCharges[0]; 
                    return true; 
                }
                if(validCharges.Count > 1)
                {
                    // if there are more than 1 zeroes found,
                    // the correct charge state will be the closest to the tested peak
                    nextPeak = validCharges.MinBy(i => i.Item2 - maxPeak.Mz);
                    return true; 
                }
            }
            // returns if the count is less than 1. 
            nextPeak = null;
            return false; 
        }
        public void FindAllRelatedChargeStates((int, double) chargeMzPair, List<MzPeak> peaks, 
            double peakFindingThreshold, List<(int, double)> chargeMzPairs)
        {
            // two while loops: one for grabbing if higher, one for grabbing if lower
            double minMzVal = peaks.MinBy(i => i.Mz).Mz;
            double maxMzVal = peaks.MaxBy(i => i.Mz).Mz;

            double currentMz = chargeMzPair.Item2;

            // initial charge state + 1 is the current for the first iteration
            int firstFoundChargeState = chargeMzPair.Item1;
            int currentCharge = firstFoundChargeState + 1;

            // runs while loop for lower m/z values
            FindAllChargeStateEnvelopePeaks(currentMz, currentCharge, minMzVal, peakFindingThreshold,
                peaks, chargeMzPairs, CalculatePlusOneZValue);
            // runs while loop for higher m/z values
            FindAllChargeStateEnvelopePeaks(currentMz, currentCharge, minMzVal, peakFindingThreshold,
                peaks, chargeMzPairs, CalculateMinusOneZValue, incrementZ: false);

            
        }
        /// <summary>
        /// Runs a while loop that identifies related charge states. Either decrements or increments 
        /// the currentZ value. 
        /// </summary>
        /// <param name="currentMz"></param>
        /// <param name="currentZ"></param>
        /// <param name="lastMZ"></param>
        /// <param name="peakFindingThreshold"></param>
        /// <param name="peaks"></param>
        /// <param name="chargeMzPairs"></param>
        /// <param name="nextValueFunction"></param><summary>Function used to calculate the next m/z value.</summary>
        /// <param name="incrementZ"></param><summary>True is accessing values below the first m/z. False gets above first m/z.</summary>
        public void FindAllChargeStateEnvelopePeaks(double currentMz, int currentZ, 
            double lastMZ, double peakFindingThreshold, List<MzPeak> peaks, List<(int, double)> chargeMzPairs, 
            Func<double, int, double> nextValueFunction, bool incrementZ = true)
        {
            int missedChargeState = 0;
            double previousMz = 0D;
            // crossedLastMz only evaluates to true when the last values are exceeded, so need to evaluate the opposite 
            // in order for the loop to continue. 
            int iterationSafegaurd = 0; 
            while (!CrossedLastMz(currentMz, previousMz, lastMZ))
            {
                // calculate next Mz

                    // filter
                    CalculateThresholds(currentMz, peakFindingThreshold, out double minVal, out double maxVal);
                IEnumerable<MzPeak> filteredPeaks = peaks.Where(peak => peak.Mz >= minVal && peak.Mz <= maxVal).ToList();
                // check if filtered count is greater than zero
                switch (filteredPeaks.Count())
                {
                    // count 0 means there were no matched charge states, but we also want to be okay with 
                    // missing a single charge state, so we increment a count that will break the loop if it goes above 
                    // one. 
                    case 0: missedChargeState++; break;
                    // Adds the currentZ and the filtered peaks if the matched peaks count is 1. 
                    case 1: chargeMzPairs.Add((currentZ, filteredPeaks.First().Mz)); break;
                    // if count (more than one matched peak) is greater than 1, we get the closest to the theoretical m/z value. 
                    case > 1: chargeMzPairs.Add((currentZ, filteredPeaks.MinBy(i => i.Mz - currentMz).Mz)); break;
                }
                // break the entire method if we go above one missed charge state.
                if (missedChargeState > 1) return;

                previousMz = currentMz; 
                currentMz = nextValueFunction(previousMz, currentZ);
                if (incrementZ) currentZ++;
                else currentZ--;

                iterationSafegaurd++;

                // break in the event that I accidentally spawned an infinite loop. 
                if (iterationSafegaurd > peaks.Count())
                {
                    throw new ArgumentException("Accidentally generated infinite loop.");
                }
            }
        }
        // relatively complicated logical construct that is used by the while loop to determine continuation. 
        // Checks the current and previous values against the last specified value. If the last value is higher than 
        // the initial value, the current and previous values will be greater and less than the last value, respectively. 
        // If the last value is lower than the initial value, the current and previous values will be less than and greater than 
        // the last value, respectively. Returns false while both of these are false, but returns true when either one condition or the other is met. 
        // It is logically impossible for both sides to be true at the same time. So to continue, while loop needs to evaluate to the opposite. 
        public bool CrossedLastMz(double current, double previous, double last)
        {
            return ((current < last) && (previous > last)) || ((current > last) && (last > previous)); 
        }
        public void CalculateThresholds(double value, double tolerance, out double minVal, out double maxVal)
        {
            minVal = value - tolerance; 
            maxVal = value + tolerance;
        }
        public bool IsWithinThreshold(double value, double lowVal, double highVal)
        {
            return value > lowVal && value < highVal; 
        }
        public double CalculatePlusOneZValue(double currentMZ, int currentZ)
        {
            return currentMZ - Math.Log(currentZ + 1); 
        }
        public double CalculateMinusOneZValue(double currentMZ, int currentZ)
        {
            return currentMZ + Math.Log(currentZ - 1); 
        }
        public void UpdatePeakList(List<MzPeak> peaksList, List<double> mzValsToRemove)
        {
            peaksList.RemoveAll(item => mzValsToRemove.Contains(item.Mz)); 
        }
        public bool IsConsecutiveChargeState(double peakMz1, double peakMz2, 
            int charge, double zeroThreshold)
        {
            double result = ChargeStateCalculation(peakMz1, peakMz2, charge);
            double maxTol = 0D + zeroThreshold;
            double minTol = 0D - zeroThreshold; 
            if(result > minTol && result < maxTol)
            {
                return true;
            }
            else
            {
                return false; 
            }
        }
        public double ChargeStateCalculation(double peakMz1, double peakMz2, 
            int charge)
        {
            // if the peaks are in the same charge state envelope, then multiplying two peaks 
            // by their two correct integers, then their difference will be zero: 
            // m/z * z = m/(z+1) * (z+1)
            return peakMz1 * charge - (peakMz2 * charge + 1); 
        }
        private void CreatesChargesAndLogZ(int minZ, int maxZ)
        {
            int numberOfCharges = maxZ - minZ;
            Charges = new int[numberOfCharges + 1]; 
            LogZ= new double[numberOfCharges + 1];
            for(int i = 0; i < numberOfCharges; i++)
            {
                Charges[i] = minZ + i;
                LogZ[i] = Math.Log((double)(minZ + i)); 
            }
        }
    }
    public class ChargeStateEnvelope
    {
        public List<ChargeState> ChargeStates { get; set; }
        private List<MzPeak> InitialPeaks { get; set; }
        public ChargeStateEnvelope()
        {
            InitialPeaks = new List<MzPeak>();
            ChargeStates = new List<ChargeState>(); 
        }
    }
    public class ChargeState
    {

    }
}
