using Thermo.TNG.Factory;
using Thermo.Interfaces.FusionAccess_V1;
using Thermo.Interfaces.FusionAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.Control;
using Thermo.Interfaces.InstrumentAccess_V1.Control.Acquisition;


// Wrapper for the thermo provided API. 
namespace InstrumentControl
{
	// Implements the wrapper around the Thermo API. 
	public class FusionLumosAPI
	{
		public string InstrumentID { get; set; }
		public string InstrumentName { get; set; }

		// Each of the below properties is a different facet of the instrument API. 
		public IFusionInstrumentAccessContainer InstAccessContainer { get; private set; }
		public IFusionInstrumentAccess InstAccess { get; private set; }
		public IFusionMsScanContainer MSScanContainer { get; private set; }
		public IAcquisition InstAcq { get; private set; }
		public IControl InstControl { get; private set; }

		// Constructor for the API. Initially only sets the InstAccessContainer property.
		// Other properties are filled when calling GetInstAccess(). This is because
		// Factory<IFusionInstrumentAccessContainer>.Create() starts online service,
		// but doesn't start the instrument connection. GetInstAccess() actually starts the
		// instrument connection. 
		public FusionLumosAPI()
		{
			InstAccessContainer = Factory<IFusionInstrumentAccessContainer>.Create(); 
		}
		// values depend on event handling.
		// TODO: connect ServiceConnected and InstrumentConnected to their
		// corresponding component of the instrument of the Thermo API. 
		public bool ServiceConnected { get; private set; }
		public bool InstrumentConnected { get; private set; }

		// closes the instrument connection.
		// marked internal 
		internal void StartOnlineAccess()
		{
			// TODO: Add error handling if connection can't be reached. 
			// Use a try catch, return error messaging. 
			InstAccessContainer.StartOnlineAccess(); 
		}
		// Closes the instrument connection. 
		internal void CloseConnection()
		{
			InstAccessContainer.Dispose(); 
		}
		// GetInstAccess beings the instrument connection and fills the other properties.
		internal void GetInstAccess(int p) // unsure what int p is for. API documentation doesn't have any info on it
		{
			InstAccess = InstAccessContainer.Get(p);
			// do not change order. InstAccess must be filled first as the other
			// properties depend on it to be filled themselves.
			InstControl = InstAccess.Control;
			InstAcq = InstControl.Acquisition;
			InstrumentID = InstAccess.InstrumentId.ToString();
			InstrumentName = InstAccess.InstrumentName; 
		}
	}
}

