using System;
using System.Collections.Generic;
using System.Text;
using Proteomics;
using UsefulProteomicsDatabases;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Text.RegularExpressions;

namespace Search_pool
{
    public class SearchPool
    {
        public SearchPool()
        {

        }

        public Parameters Parameters { get; set; }

        public static Dictionary<char, char> DegenerateTable = new Dictionary<char, char>
        {
            {'G', 'G'},{'A', 'G'},{'S', 'G'},{'C', 'G'},{'V', 'G'},{'T', 'G'},{'I', 'G'},
            {'g', 'G'},{'a', 'G'},{'s', 'G'},{'c', 'G'},{'v', 'G'},{'t', 'G'},{'i', 'G'},
            {'P', 'P'},
            {'p', 'P'},
            {'L', 'L'},{'D', 'L'},{'N', 'L'},{'E', 'L'},{'Q', 'L'},{'M', 'L'},
            {'l', 'L'},{'d', 'L'},{'n', 'L'},{'e', 'L'},{'q', 'L'},{'m', 'L'},
            {'K', 'K'},{'R', 'K'},
            {'k', 'K'},{'r', 'K'},
            {'H', 'Y'},{'F', 'Y'},{'Y', 'Y'},
            {'h', 'Y'},{'f', 'Y'},{'y', 'Y'},
            {'W', 'W'},
            {'w', 'W'},
            {'X', 'X'},
            {'x', 'X'},

            {'U', 'G'},
            {'u', 'G'},
            {'B', 'L'},{'Z', 'L'},
            {'b', 'L'},{'z', 'L'},
        };

        public static Dictionary<string, char> ResiduesDictionary = new Dictionary<string, char>
        {
            {"ALA", 'A'},
            {"ARG", 'R'},
            {"ASN", 'N'},
            {"ASP", 'D'},
            {"CYS", 'C'},
            {"GLU", 'E'},
            {"GLN", 'Q'},
            {"GLY", 'G'},
            {"HIS", 'H'},
            {"ILE", 'I'},
            {"LEU", 'L'},
            {"LYS", 'K'},
            {"MET", 'M'},
            {"PHE", 'F'},
            {"PRO", 'P'},
            {"PYL", 'O'},
            {"SEL", 'U'},
            {"SER", 'S'},
            {"THR", 'T'},
            {"TRP", 'W'},
            {"TYR", 'Y'},
            {"VAL", 'V'},
            {"MSE", 'X'}  //print_seq = ['phenix.print_sequence', params.queryfile, '--letter_for_mse=X']
        };

        public static List<Protein> LoadProteinDb(string fileName, bool generateTargets, DecoyType decoyType, List<string> localizeableModificationTypes, bool isContaminant, int MaxThreadsToUsePerFile)
        {
            List<string> dbErrors = new List<string>();
            List<Protein> proteinList = new List<Protein>();

            string theExtension = Path.GetExtension(fileName).ToLowerInvariant();
            bool compressed = theExtension.EndsWith("gz"); // allows for .bgz and .tgz, too which are used on occasion
            theExtension = compressed ? Path.GetExtension(Path.GetFileNameWithoutExtension(fileName)).ToLowerInvariant() : theExtension;

            if (theExtension.Equals(".fasta") || theExtension.Equals(".fa"))
            {
                proteinList = ProteinDbLoader.LoadProteinFasta(fileName, generateTargets, decoyType, isContaminant, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, out dbErrors, MaxThreadsToUsePerFile);
            }

            return proteinList.Where(p => p.BaseSequence.Length > 0).ToList();
        }

        public static string DegenerateSequence(string sequence)
        {
            string dSequence = "";
            foreach (var c in sequence)
            {
                dSequence += DegenerateTable[c];
            }
            return dSequence;
        }
        public static Protein DegenerateProtein(Protein protein)
        {
            var updateBaseSeq = DegenerateSequence(protein.BaseSequence);

            Protein dProtein = new Protein(protein, updateBaseSeq);
            return dProtein;
        }
        public static Protein[] DegenerateProteinDa(List<Protein> proteins, int MaxThreadsToUsePerFile)
        {
            Protein[] dProteins = new Protein[proteins.Count];

            int[] threads = Enumerable.Range(0, MaxThreadsToUsePerFile).ToArray();

            Parallel.ForEach(threads, (index) =>
            {
                for (; index < proteins.Count; index += MaxThreadsToUsePerFile)
                {
                    dProteins[index] = DegenerateProtein(proteins[index]);
                }
            });

            return dProteins;
        }

        public static List<string> LoadQueryMap(string fileName)
        {
            List<string> queries = new List<string>();

            string currentSeq = "";

            HashSet<string> seenAminoAcids = new HashSet<string>();

            HashSet<string> SeenSeqId = new HashSet<string>();

            if (Path.GetExtension(fileName) == ".pdb")
            {
                using (StreamReader streamReader = new StreamReader(fileName))
                {
                    while (streamReader.Peek() != -1)
                    {
                        string line = streamReader.ReadLine();

                        if (line.StartsWith("TER") || line.StartsWith("END"))
                        {
                            queries.Add(currentSeq.ToString());
                            currentSeq = "";
                        }
                      
                        if (line.StartsWith("ATOM") || line.StartsWith("HETATM"))
                        {                           
                            var sp = Regex.Split(line, @"\s+");
                            if (SeenSeqId.Contains(sp[4]))
                            {
                                string key = sp[4] + "_" + sp[5];
                                if (!seenAminoAcids.Contains(key))
                                {
                                    seenAminoAcids.Add(key);                                  
                                    currentSeq += ResiduesDictionary[sp[3]];
                                }
                            }
                            else
                            {                         
                                if (SeenSeqId.Count > 0)
                                {
                                    queries.Add(currentSeq.ToString());
                                    currentSeq = "";                          
                                }
                                else
                                {
                                    string key = sp[4] + "_" + sp[5];
                                    if (!seenAminoAcids.Contains(key))
                                    {
                                        seenAminoAcids.Add(key);
                                        currentSeq += ResiduesDictionary[sp[3]];
                                    }
                                }
                                SeenSeqId.Add(sp[4]);
                            }

                        }
                    }
                }
            }

            return queries;
        }

        public static List<string> DegenerateQueryMap(List<string> queries)
        {
            List<string> dQueries = new List<string>();

            foreach (var q in queries)
            {
                dQueries.Add( DegenerateSequence(q));
            }

            return dQueries;
        }
   
        
    }
}
