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
            Parameters = new Parameters();
        }

        public Parameters Parameters { get; private set; }

        public void Run(string db_path, string query_path, int matrixId)
        {
            var queries = DataBaseManipulation.LoadQueryMap(query_path);

            var dQueries = DataBaseManipulation.DegenerateQueryMap(queries);

            var proteins = DataBaseManipulation.LoadProteinDb(db_path, true, DecoyType.None, Parameters.MaxThreadsToUse);

            var dProteins = DataBaseManipulation.DegenerateProteinDa(proteins, Parameters.MaxThreadsToUse);

            int[][] scores = new int[dQueries.Count][];

            string[][] matchedSequences = new string[dQueries.Count][];

            int gap = 32767;

            for (int i = 0; i < dQueries.Count; i++)
            {
                scores[i] = new int[dProteins.Length];

                matchedSequences[i] = new string[dProteins.Length];

                int[] threads = Enumerable.Range(0, Parameters.MaxThreadsToUse).ToArray();                

                Parallel.ForEach(threads, (index) =>
                {
                    for (; index < dProteins.Length; index += Parameters.MaxThreadsToUse)
                    {
                        MatrixAlignment ma = new MatrixAlignment(dProteins[index].BaseSequence, dQueries[i]);

                        if (matrixId == 0)
                        {
                            MatrixAlignment.AlignMatrix(ma, gap, MatrixAlignment.ScoreTable);
                        }
                        else if (matrixId == 1)
                        {
                            MatrixAlignment.AlignMatrix_UpLeft(ma, gap, MatrixAlignment.ScoreTable);
                        }       

                        scores[i][index] = ma.MaxScore;

                        var path = MatrixAlignment.Traceback(ma);

                        var seqCom = MatrixAlignment.SequenceComparison(ma, path);

                        matchedSequences[i][index] = string.Join('-', seqCom.Select(p => new string(p)));
                    }
                });
            }

            int[] addup_score = new int[dProteins.Length];
            for (int j = 0; j < dProteins.Length; j++)
            {
                for (int i = 0; i < dQueries.Count; i++)
                {
                    addup_score[j] += scores[i][j];
                }
            }

            int[] indexes = Enumerable.Range(0, dProteins.Length).ToArray();
            Array.Sort(addup_score, indexes);
            Array.Reverse(addup_score);
            Array.Reverse(indexes);

            List<Target> targets = new List<Target>();

            int ind_add = 0;
            foreach (var ind in indexes.Take(20))
            {
                var target = Target.SetTarget(addup_score[ind_add], dProteins[ind].Accession);
                ind_add++;
                target.SingleScores = new int[dQueries.Count];
                target.Matched_sequences = new string[dQueries.Count];
                for (int i = 0; i < dQueries.Count; i++)
                {
                    target.SingleScores[i] = scores[i][ind];
                    target.Matched_sequences[i] = matchedSequences[i][ind];
                }
                targets.Add(target);
            }

            //Write candidate out.
            var fileName = "target_matrix" + matrixId + ".tsv";
            string outFilePath = Path.Combine(Path.GetDirectoryName(query_path), fileName);
            WriteTargets(targets, outFilePath);
        }

        public static void WriteTargets(List<Target> targets, string outFilePath)
        {
            using (StreamWriter writer = new StreamWriter(outFilePath))
            {
                foreach (var t in targets)
                {
                    writer.WriteLine(t.Protein + "\t" + t.Addup_Score + "\t");
                    for (int i = 0; i < t.SingleScores.Length; i++)
                    {
                        writer.WriteLine("\t" + t.SingleScores[i] + "\t" + t.Matched_sequences[i]);
                    }
                }
            }
        }
    }
}
