using System;
using System.Collections;
using System.Linq;
using MathNet;
using System.Collections.Generic;
using System.Threading.Tasks;
using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;
using OxyPlot.Reporting;


namespace FtmsSimulator
{
	public class InstrumentParameters
	{
		public bool Decaying { get; set; }
		public double ChargeConstant { get; set; }
		public double DeltaZ { get; set; }
		public double Lambda { get; set; }
		public double FinalIonPercentage { get; set; }

		public InstrumentParameters(double deltaZ, double lambda, double chargeConstant, bool decaying, double finalIonPercentage)
		{
			DeltaZ = deltaZ;
			Lambda = lambda;
			ChargeConstant = chargeConstant; 
			Decaying = decaying;
			FinalIonPercentage = finalIonPercentage; 
		}
	}
	public class Transient
	{
		double MinT { get; set; }
		double MaxT { get; set; }
		double DeltaT { get; set; }
		public double[] SimulatedValues { get; set; }
		int Q { get; set; }
		double M { get; set; }
		bool Decaying { get; set; }
		double DecayConstant { get; set; }
		int NumberOfIons { get; set; }
		public Transient(double mass, int charge, double chargeConstant, 
			int numberOfIons, double deltaZ, double lambda, double minT, 
			double maxT, double deltaT, bool decaying, bool providedDecayConstant, double finalIonPercentage, 
			double? decayConstant = null)
		{
			M = mass;
			Q = charge;
			Decaying = decaying;
			MinT = minT;
			MaxT = maxT;
			DeltaT = deltaT;
			Decaying = decaying;
			NumberOfIons = numberOfIons; 
			if (Decaying) 
			{
				if (providedDecayConstant)
				{
					DecayConstant = decayConstant.Value;
				}
				else
				{
					DecayConstant = CalculateIonDecayConstant(maxT: maxT, percentDecay: finalIonPercentage, initialIons: numberOfIons);
				}				
			}
			else
			{
					DecayConstant = 0.0;
			}
			SimulatedValues = SimulateTransient(M, Q, chargeConstant, numberOfIons, deltaZ, lambda, MinT, MaxT, DeltaT, DecayConstant); 
		}
		public Transient((int, double, double) ionCountTuple, InstrumentParameters instrPara, 
			double minT, double maxT, double deltaT, bool providedDecayConstant, double finalIonPercentage,
			double? decayConstant = null)
		{
			M = ionCountTuple.Item2;
			Q = ionCountTuple.Item1;
			Decaying = instrPara.Decaying;
			MinT = minT;
			MaxT = maxT;
			DeltaT = deltaT;
			NumberOfIons = (int)Math.Round(ionCountTuple.Item3);

			if (Decaying)
			{
				if (providedDecayConstant)
				{
					DecayConstant = decayConstant.Value;
				}
				else
				{
					DecayConstant = CalculateIonDecayConstant(maxT: maxT, percentDecay: finalIonPercentage, initialIons: NumberOfIons);
				}
			}
			else
			{
				DecayConstant = 0.0;
			}
			SimulatedValues = SimulateTransient(M, Q, instrPara.ChargeConstant, NumberOfIons, 
				instrPara.DeltaZ, instrPara.Lambda, MinT, MaxT, DeltaT, DecayConstant);
		}

		public static double[] SimulateTransient(double mass, int chargeInt, double chargeConstant, int numberOfIons, double deltaZ, double lambda, double minT, double maxT, double deltaT, double decayConstant)
		{
			// need to convert mass from Da to kg: 
			double massKg = mass * 1.66054E-27;
			double charge = (double)chargeInt * 1.602E-19; 
			// I(t) = -q*N*e^(-at)*w*dz/lambda*sin(w*t)
			// Calculate time array
			double[] timeArray = CalculateTimeArray(minT, maxT, deltaT);
			double[] transientResult = new double[timeArray.Length]; 
			// Calculate omega: 
			
			double omega = Math.Sqrt((charge / massKg) * chargeConstant);
			// Create an array containing the number of ions at each timepoint t. 
			int[] ionCountArray = CalculateInstantaneousIonCount(decayConstant, timeArray, numberOfIons); 

			// iterate through the time series array to calculate image current at each timepoint. 
			for(int i = 0; i < timeArray.Length; i++)
			{
				transientResult[i] = CalculateInstantaneousTransient(charge, (double)ionCountArray[i], deltaZ, lambda, omega, timeArray[i]); 
			}
			return transientResult; 
		}

		public static double CalculateInstantaneousTransient(double charge, double ionCount, double deltaZ, double lambda, double omega, double time)
		{
			double instantaneousCharge = -charge * ionCount * omega * deltaZ / lambda * Math.Sin(omega * time);
			return instantaneousCharge; 
		}
		
		public static double[] CalculateTimeArray(double minT, double maxT, double deltaT)
		{
			double totalTimeSteps = (maxT - minT) / deltaT;
			double[] timeArray = Enumerable.Range(0, (int)totalTimeSteps).Select(i => (double)i * deltaT).ToArray();
			return (timeArray);
		}
		public static double CalculateIonDecayConstant(double maxT, double percentDecay, int initialIons)
		{
			int detectionEndIonCount = (int)(initialIons * percentDecay);
			double decayConstant = Math.Log(initialIons) / maxT - Math.Log(detectionEndIonCount) / maxT;
			return decayConstant; 
		}
		public static int[] CalculateInstantaneousIonCount(double ionDecayConstant, double[] timeArray, double initialIons)
		{
			int[] instnIonCount = new int[timeArray.Length]; 
			for(int i = 0; i < timeArray.Length; i++)
			{
				instnIonCount[i] = (int)(initialIons * Math.Pow(Math.E, -ionDecayConstant * timeArray[i])); 
			}
			return instnIonCount; 
		}
		public static List<(int, double, double)> GenerateIonCounts(int[] chargeStates, double[] isotopicIntensities, double[] ionCounts, double[] masses)
		{
			if (chargeStates.Length != ionCounts.Length) throw new InvalidOperationException("chargeStates.Length must be equal to ionCounts.Length"); 
			
			List<(int, double, double)> results = new List<(int, double, double)>(); 

			// for each ion count and for each isotopic masses, multiply ion count * calculated isotopic intensity
			// 
			for(int i = 0; i < ionCounts.Length; i++)
			{
				for(int j = 0; j < masses.Length; j++)
				{
					double tempIonCount = isotopicIntensities[j] * ionCounts[i];
					results.Add((chargeStates[i], masses[j], tempIonCount)); 

				}
			}
			return results; 
		}
		public static double[] SumMultipleTransients(List<Transient> transientList)
		{
			double[] summedTransient = new double[transientList[0].SimulatedValues.Length];
			
			for(int i = 0; i < summedTransient.Length; i++)
			{
				for(int j = 0; j < transientList.Count; j++)
				{
					summedTransient[i] += transientList[j].SimulatedValues[i]; 
				}
			}
			return summedTransient; 
		}
		public static List<Transient> GenerateManyTransients(List<(int, double, double)> ionCountsList, 
			InstrumentParameters instrPara, 
			double minT, double maxT, double deltaT, 
			double ionDecayConstant)
		{
			List<Transient> transientList = new List<Transient>();
			Parallel.ForEach(ionCountsList, i => 
			{
				var tempTransient = new Transient(i, instrPara, minT, maxT, deltaT, true, 0.60, ionDecayConstant);
				transientList.Add(tempTransient); 
			});
			return transientList; 
		}
	}
	public class TransientPlot
	{
		public static PlotModel PlotTransient(double[] summedTransient)
		{
			double[] xvals = Enumerable.Range(0, summedTransient.Length).Select(i => (double)i).ToArray();
			double[] yvals = summedTransient;

			var model = new PlotModel
			{
				Title = "Transient plot"
			};
			model.Axes.Add(new LinearAxis
			{
				Position = AxisPosition.Left, 
				Title = "Current"
			}); 
			model.Axes.Add(new LinearAxis
			{
				Position = AxisPosition.Bottom, 
				Title = "Time (ms)"
			});
			var ls = new LineSeries
			{
				Title = "Transient"
			}; 
			for(int i = 0; i < xvals.Length; i++)
			{
				ls.Points.Add(new DataPoint(xvals[i], yvals[i])); 
			}
			model.Series.Add(ls);
			return model; 
		}
	}
}

