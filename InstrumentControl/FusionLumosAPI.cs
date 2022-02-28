using Thermo.TNG.Factory;
using Thermo.Interfaces.FusionAccess_V1;
using Thermo.Interfaces.FusionAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.Control;
using Thermo.Interfaces.InstrumentAccess_V1.Control.Acquisition;


// Wrapper for the thermo provided API. 
namespace InstrumentControl
{
	enum Module { None = -1, Default = 0}
	public class FusionLumosAPI
	{
		public string InstrumentID { get; set; }
		public string InstrumentName { get; set; }

		public IFusionInstrumentAccessContainer InstAccessContainer { get; private set; }
		public IFusionInstrumentAccess InstAccess { get; private set; }
		public IFusionMsScanContainer MSScanContainer { get; private set; }
		public IAcquisition InstAcq { get; private set; }
		public IControl InstControl { get; private set; }
		public FusionLumosAPI()
		{
			InstAccessContainer = Factory<IFusionInstrumentAccessContainer>.Create(); 
		}
		// values depend on event handling 
		public bool ServiceConnected { get; private set; }
		public bool InstrumentConnected { get; private set; }
		internal void StartOnlineAccess()
		{
			// TODO: Add error handling if connection can't be reached. 
			// Use a try catch, return error messaging. 
			InstAccessContainer.StartOnlineAccess(); 
		}
		internal void CloseConnection()
		{
			InstAccessContainer.Dispose(); 
		}
		internal void GetInstAccess(int p) // unsure what int p is for. API documentation doesn't have any info on it
		{
			InstAccess = InstAccessContainer.Get(1);
			InstControl = InstAccess.Control;
			InstAcq = InstControl.Acquisition;
			InstrumentID = InstAccess.InstrumentId.ToString();
			InstrumentName = InstAccess.InstrumentName; 
		}
	}
}

