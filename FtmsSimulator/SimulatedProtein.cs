using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Chemistry;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot.Reporting; 


namespace FtmsSimulator
{
	public class SimulatedProtein : AminoAcidPolymer
	{
		SimulatedProtein()
		{

		}
		public SimulatedProtein(string sequence) : base(sequence)
		{
		}

		public SimulatedProtein(AminoAcidPolymer aminoAcidPolymer)
							: this(aminoAcidPolymer, true)
		{
		}
        public SimulatedProtein(AminoAcidPolymer aminoAcidPolymer, bool includeModifications)
            : base(aminoAcidPolymer, includeModifications)
        {
            Parent = aminoAcidPolymer;
            StartResidue = 0;
            EndResidue = Length - 1;
        }

        public SimulatedProtein(AminoAcidPolymer aminoAcidPolymer, int firstResidue, int length)
                    : this(aminoAcidPolymer, firstResidue, length, true)
        {
        }

        public SimulatedProtein(AminoAcidPolymer aminoAcidPolymer, int firstResidue, int length, bool includeModifications)
                    : base(aminoAcidPolymer, firstResidue, length, includeModifications)
        {
            Parent = aminoAcidPolymer;
            StartResidue = firstResidue;
            EndResidue = firstResidue + length - 1;
            PreviousResidue = aminoAcidPolymer.GetResidue(StartResidue - 1);
            NextResidue = aminoAcidPolymer.GetResidue(EndResidue + 1);
        }
        public int StartResidue { get; set; }

        /// <summary>
        /// The amino acid number this peptide is located in its parent
        /// </summary>
        public int EndResidue { get; set; }

        /// <summary>
        /// The amino acid polymer this peptide came from. Could be null
        /// </summary>
        public AminoAcidPolymer Parent { get; set; }

        /// <summary>
        /// The preceding amino acid in its parent
        /// </summary>
        public Residue PreviousResidue { get; set; }

        /// <summary>
        /// The next amino acid in its parent
        /// </summary>
        public Residue NextResidue { get; set; }

        public IEnumerable<SimulatedProtein> GenerateAllModificationCombinations()
        {
            // Get all the modifications that are isotopologues
            var isotopologues = GetUniqueModifications<ModificationWithMultiplePossibilitiesCollection>().ToArray();

            // Base condition, no more isotopologues to make, so just return
            if (isotopologues.Length < 1)
            {
                yield break;
            }

            // Grab the the first isotopologue
            ModificationWithMultiplePossibilitiesCollection isotopologue = isotopologues[0];

            // Loop over each modification in the isotopologue
            foreach (OldSchoolModification mod in isotopologue)
            {
                // Create a clone of the peptide, cloning modifications as well.
                SimulatedProtein simProt = new SimulatedProtein(this);

                // Replace the base isotopologue mod with the specific version
                simProt.ReplaceModification(isotopologue, mod);

                // There were more than one isotopologue, so we must go deeper
                if (isotopologues.Length > 1)
                {
                    // Call the same rotuine on the newly generate peptide that has one less isotopologue
                    foreach (var subpeptide in simProt.GenerateAllModificationCombinations())
                    {
                        yield return subpeptide;
                    }
                }
                else
                {
                    // Return this peptide
                    yield return simProt;
                }
            }
        }

        public SimulatedProtein GetSubProt(int firstResidue, int length)
        {
            return new SimulatedProtein(this, firstResidue, length);
        }

        public void CalculateIsotopicDistribution()
		{
            ChemicalFormula formula = GetChemicalFormula();

		}
		public void CreateChargeStateProbabilityDistribution(string distributionType,
			double mass,
			int minCharge,
			int maxCharge)
		{
			// Calculate binomial and negative binomial distributions

			//Binomical distribution: 
			// 
		}

        public static int Factorial(int x)
		{
            // max value before it returns NaN is 2,147,483,647
            if (x <= 1)
			{
                return 1;
			}
			else
			{
                return x * Factorial(x - 1); 
			}
		}
        public static double BinomialProbability(int trials, int successes, double probabilityOfSuccess)
		{
            int frontTerm1 = Factorial(trials); 
            int frontTerm2 = Factorial((trials - successes)) * Factorial(successes);
            double frontTerm = (double)frontTerm1 / (double)frontTerm2; 
            double successProb = Math.Pow(probabilityOfSuccess, (double)successes); 
            double unsuccessProb = Math.Pow((1 - probabilityOfSuccess), (double)(trials - successes));
            double result = frontTerm * successProb * unsuccessProb;
            return result; 
		}
        public static double[] CalculateBinomialDistribution(int[] seqOfSuccesses, 
            int trials, 
            double probabilityOfSuccesses)
		{
            double[] resultsArray = new double[seqOfSuccesses.Length]; 
            for(int i = 0; i < seqOfSuccesses.Length; i++)
			{
                resultsArray[i] = BinomialProbability(trials, seqOfSuccesses[i], probabilityOfSuccesses); 
			}
            return resultsArray; 
		}
        public static double[] CreateLognormalDistribution(double shape, double scale, double location, double[] valueArray)
        {
            // for lognormal formula information and equation: https://www.statisticshowto.com/lognormal-distribution/
            double[] resultsArray = new double[valueArray.Length];
            double frontTerm = 1 / (shape * Math.Sqrt(2 * Math.PI));
            for (int i = 0; i < valueArray.Length; i++)
			{
                double expTerm = -(Math.Pow(Math.Log(i - location), 2) / (2 * Math.Pow(shape, 2)));
                resultsArray[i] = frontTerm * expTerm; 
            }
            return resultsArray; 
        }
        public static double[] CreateWeibullDistribution(double shape, double scale, double location, double[] valueArray)
		{
            double[] resultsArray = new double[valueArray.Length];
            double frontTerm = shape / scale; 

            for(int i = 0; i < valueArray.Length; i++)
			{
                double term1 = (valueArray[i] - location) / scale;
                double term2 = -(Math.Pow(term1, shape));
                resultsArray[i] = frontTerm * Math.Pow(term1, (shape - 1)) * Math.Exp(term2); 
			}
            return resultsArray; 

		}
	}
	public class SimulatedProteinPlot
	{
        public static PlotModel PlotIsotopicDistribution(IsotopicDistribution isoEnvelope)
		{
            double[] xvals = isoEnvelope.Masses.ToArray();
            double[] yvals = isoEnvelope.Intensities.ToArray();
            var model = new PlotModel
            {
                Title = "Isotopic Distribution Plot"
            };
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left
            });
            model.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Bottom
            });
            var ls = new LineSeries
            {
                Title = "Distribution"
            }; 
            for(int i = 0; i < xvals.Length; i++)
			{
                ls.Points.Add(new DataPoint(xvals[i], yvals[i]*10E9)); 
			}
            model.Series.Add(ls);
            return model; 

		}
	}
}
