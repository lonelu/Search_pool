using NUnit.Framework;
using Search_pool;
using System.IO;
using System.Linq;
using System;


namespace Test
{
    public class UnitTest_EValue
    {
        [Test]
        public void Test_SpougeStoE()
        {
            var kbp = EValue.Set_Blast_KarlinBlk(0);

            var gbp = EValue.Set_Blast_GumbelBlk(32767, 32767, 0);

            var e = EValue.BLAST_SpougeStoE(64, kbp, gbp, 20, 20);

            var e1 = EValue.BLAST_KarlinStoE_simple(64, kbp, 20);
        }
    }
}
