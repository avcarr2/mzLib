using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.ReaderFactories
{
    internal class BrukerReaderFactory : BaseReaderFactory, IReaderFactory
    {
        public MsDataFile Reader { get; private set; }

        internal BrukerReaderFactory(string filePath) : base(filePath)
        {
            FilePath = filePath;
            Reader = CreateReader();
        }

        public MsDataFile CreateReader()
        {
            return new Bruker.Bruker(FilePath); 
        }
    }
}
