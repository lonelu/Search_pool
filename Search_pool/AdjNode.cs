using System;
using System.Collections.Generic;
using System.Text;

namespace Search_pool
{
    public struct AdjNode
    {
        public int MaxCost { get; set; }

        //Source Tuple<int, int, int>(pos j, pos i, from 0:lu|1:l|-1:u )
        public Tuple<int, int, int> Source { get; set; }
    }
}
