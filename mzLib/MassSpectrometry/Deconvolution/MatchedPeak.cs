using System;
using Chemistry;

namespace MassSpectrometry;

public class MatchedPeak
{
    public double Mz { get; set; }
    public double Intensity { get; set; }
    public double Charge { get; set; }
    public double MatchError { get; set; }
    public double NeutralMass => CalculateNeutralMass();
    public MatchedPeak(double mz, double intensity, double charge)
    {
        Mz = mz;
        Intensity = intensity;
        Charge = charge;
    }

    public void CalculateMatchError(ChargeStateLadder theoreticalLadder)
    {
        int index = ChargeStateIdentifier.GetBucket(theoreticalLadder.MzVals, Mz);
        MatchError = Math.Abs((Mz - theoreticalLadder.MzVals[index]) / Mz * 1E6);
    }

    private double CalculateNeutralMass()
    {
        return Mz.ToMass((int)Math.Round(Charge));
    }
}