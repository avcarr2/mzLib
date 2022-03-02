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
			RealTimeSampleProcessingExample dataReceiver = new(); 

			// main loop 
			while (api.InstAccessContainer.ServiceConnected)
			{
				// MsScanArrived is an event.
				// The += means "subscribes to".
				// So the data.Receiver.MSScanContainer_MsScanArrived method is
				// subscribed to the MsScanArrived event and will be called when the event is raised.

				/* Technically, multiple methods can subscribe to the same event, 
				 * however, it is likely doing so would require multithreading. 
				 * You would need to have multiple objects working with a single IMsScan object, 
				 * and I'm not sure how you would pass a single IMsScan to multiple threads yet
				 */
				api.MSScanContainer.MsScanArrived += dataReceiver.MSScanContainer_MsScanArrived; 
			}
			
		}
		// Instantiates the FusionLumosAPI and starts both online access and instrument access
		private static FusionLumosAPI InitializeConnection()
		{
			FusionLumosAPI api = new();
			api.StartOnlineAccess(); // call to start online connection. Need to add event and error handling. 
			api.GetInstAccess(1); // access instrument and fills the FusionLumosAPI class properties. 
			return api; 
		} 
	}
}
