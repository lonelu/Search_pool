using NUnit.Framework;
using Search_pool;
using System.IO;
using System.Linq;
using System;

namespace Test
{
    public class Tests
    {
        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void Test_LoadQueryMap()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\UnkwnPro_queries_optimized.pdb");

            var queries = SearchPool.LoadQueryMap(filePath);

            Assert.That(queries.Count == 3);

            var dQueries = SearchPool.DegenerateQueryMap(queries);

            //Assert.That(dQueries.First() == "ZGGMYZZZZGGYPZYXGYGZ");

            Assert.That(dQueries.First() == "LGGKYLLLLGGYPLYXGYGL");
        }

        [Test]
        public void Test_AlignMatrix()
        {
            string Q8I2J3 = "MDKKAREYAQDALKFIQRSGSNFLACKNLKERLENNGFINLSEGETWNLNKNEGYVLCKE" +
                "NRNICGFFVGKNFNIDTGSILISIGHIDSCALKISPNNNVIKKKIHQINVECYGSGLWHT" +
                "WFDRSLGLSGQVLYKKGNKLVEKLIQINKSVLFLPSLAIHLQNRTRYDFSVKINYENHIK" +
                "PIISTTLFNQLNKCKRNNVHHDTILTTDTKFSHKENSQNKRDDQMCHSFNDKDVSNHNLD" +
                "KNTIEHLTNQQNEEKNKHTKDNPNSKDIVEHINTDNSYPLLYLLSKELNCKEEDILDFEL" +
                "CLMDTQEPCFTGVYEEFIEGARFDNLLGSFCVFEGFIELVNSIKNHTSNENTNHTNNITN" +
                "DINDNIHNNLYISIGYDHEEIGSLSEVGARSYCTKNFIDRIISSVFKKEIHEKNLSVQEI" +
                "YGNLVNRSFILNVDMAHCSHPNYPETVQDNHQLFFHEGIAIKYNTNKNYVTSPLHASLIK" +
                "RTFELYYNKYKQQIKYQNFMVKNDTPCGSTVGSMVAANLSMPGIDIGIPQLAMHSIREIA" +
                "AVHDVFFLIKGVFAFYTYYNQVLSTCVHDK";

            string protein = SearchPool.DegenerateSequence(Q8I2J3);

            string query = "LGGKYLLLLGGYPLYXGYGL";
            int gap = 32767;

            MatrixAlignment ma = new MatrixAlignment(protein, query);

            MatrixAlignment.AlignMatrix(ma, gap, MatrixAlignment.ScoreTable);

            var path = MatrixAlignment.Traceback(ma);

            var seqCom = MatrixAlignment.SequenceComparison(ma, path);

        }


        [Test]
        public void Test_AlignMatrix2()
        {
            string O77390 = "MNKMASTHNEIIPRLGFEEMRNEMNKYGVEINQSTLKNPSTEDIQGIYSLCIKYILNKDI" +
                "QNIRIEEYTGDLKSSLPTVDGLQILPNEGKNHLQAIGNLRFLRHCEKINKILNLDNILSY" +
                "IFKPVGSHMTKLINAFIHFMKYRDQLYNENGEKIKSIQEKKNEYDVLENEYDALENELNK" +
                "LLLKHEDIRNNIINEKNIKRNYEEDIIKNQNLLNSQQSLIISLNSTKDKIVNETNELIFQ" +
                "YSRYRQKKEDLEDQIVPSPEKLQKYNEELKDHLYEHIAQFEDDRKKNEDIKNKINIADIC" +
                "IKKLVDLLTALNEHIEHTIKLHIEKKNNLQTIEKQYKSLTNEKQNFITKNTEQDKIIKET" +
                "KEFLQQEQTKWNQKIKQEQHNTILIQQKVKDIYQNVDDLNIKTNREINQINNIIKHIQDI" +
                "INHYNKNILLITELIQNTKNSHSILTHKVLNNIQKDISANM";

            string protein = SearchPool.DegenerateSequence(O77390);

            string protein_select = "YGGLGYKGLLLGLKLGGG";

            //string query = "KGKGYGLLGGXYGXKGGGPYLGGLGLKGK";
            string query = "GGGKGYGGLLYXLKLGGGGY";
            
            int gap = 32767; //Default in Blastp gapopen 

            MatrixAlignment ma = new MatrixAlignment(protein, query);

            MatrixAlignment.AlignMatrix(ma, gap, MatrixAlignment.ScoreTable);

            var path = MatrixAlignment.Traceback(ma);

            var seqCom = MatrixAlignment.SequenceComparison(ma, path);

            Assert.That(new string(seqCom[0]) == "GGKGYGGLLYXLKLGGG");
            Assert.That(new string(seqCom[1]) == "GGLGYKGLLLGLKLGGG");

            //MatrixAlignment.PrintMatrix(ma, Path.Combine(TestContext.CurrentContext.TestDirectory, @"matrix_alignment_debug.txt"));

        }

        

    }
}