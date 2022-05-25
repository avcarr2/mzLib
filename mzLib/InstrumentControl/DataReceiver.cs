using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Thermo.Interfaces.FusionAccess_V1;
using Thermo.Interfaces.FusionAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.MsScanContainer;
using Thermo.Interfaces.InstrumentAccess_V1.Control.Scans;
using Thermo.Interfaces.SpectrumFormat_V1;
using Thermo.Interfaces.FusionAccess_V1.Control.Scans;
using System.Threading;
using MassSpectrometry; 

// Abstract class used to template the individual data processing modules.
// Data processing modules should seek to free up the IMsScan yielded from e.GetScan() by processing
// it to an MsDataScan object as quickly as possible. IMsScan objects implement the IDisposable
// interface, which means they can be used with a using() statement for easier clean up. 

namespace InstrumentControl
{
	public abstract class DataReceiver
	{
		public abstract Queue<MsDataScan> ScanProcessingQueue { get; }

		public abstract void MSScanContainer_MsScanArrived(object sender, MsScanEventArgs e); 
	}
}
