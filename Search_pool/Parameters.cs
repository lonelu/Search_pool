using System;
using System.Collections.Generic;
using System.Text;

namespace Search_pool
{
    public class Parameters
    {
        public Parameters()
        {
            MaxThreadsToUsePerFile = 4;
        }

        public int MaxThreadsToUsePerFile { get; set; }
    }
}
