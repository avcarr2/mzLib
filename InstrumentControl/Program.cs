using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InstrumentControl
{
	public static class Program
	{
		public static void Main(string[] args)
		{
			FusionLumosAPI api = InitializeConnection();
			DataReceiver dataReceiver = new(); 

			while (api.InstAccessContainer.ServiceConnected)
			{
				api.MSScanContainer.MsScanArrived += dataReceiver.MSScanContainer_MsScanArrived; 
				
			}
			
		}
		public static FusionLumosAPI InitializeConnection()
		{
			FusionLumosAPI api = new();
			api.StartOnlineAccess(); // call to start online connection. Need to add event and error handling. 
			api.GetInstAccess(1); // access instrument and fills the FusionLumosAPI class properties. 
			return api; 
		} 
	}
}
