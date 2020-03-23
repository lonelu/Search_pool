using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using System.IO;

namespace Search_pool
{
    public class MatrixAlignment
    {
        public AdjNode[][] matrix { get; set; }

        public string Protein { get; set; }

        public string Query { get; set; }

        public int ProteinLength { get { return Protein.Length; } }

        public int QueryLength { get { return Query.Length; } }

        public int MaxScore { get; set; }

        public Tuple<int, int> MaxScorePos {get; set;}

        public MatrixAlignment(string protein, string query)
        {
            Protein = protein;

            Query = query;

            matrix = new AdjNode[QueryLength + 1][];
            for (int i = 0; i < QueryLength + 1; i++)
            {
                matrix[i] = new AdjNode[ProteinLength + 1];
            }
        }

        public static Dictionary<char, Dictionary<char, int>> ScoreTable = new Dictionary<char, Dictionary<char, int>>
        {
            {'G', new Dictionary<char, int>
            {
                {'G', 6 }, {'P', -6 },{'L', -5 },{'K', -8 },{'Y', -14 },{'W',-15 },{'X',-1 }
            } },
            {'P', new Dictionary<char, int>
            {
                {'G', -6 }, {'P', 8 },{'L', -4 },{'K', -8 },{'Y', -13 },{'W',-14 },{'X',-1 }
            } },
            {'L', new Dictionary<char, int>
            {
                {'G', -5 }, {'P', -4 },{'L', 6 },{'K', -5 },{'Y', -9 },{'W', -14 },{'X',-1 }
            } },
            {'K', new Dictionary<char, int>
            {
                {'G', -8 }, {'P', -8 },{'L', -5 },{'K', 11 },{'Y', -11 },{'W',-13 },{'X',-1 }
            } },
            {'Y', new Dictionary<char, int>
            {
                {'G', -14 }, {'P', -13 },{'L', -9 },{'K', -11 },{'Y', 10 },{'W',-5 },{'X',-1 }
            } },
            {'W', new Dictionary<char, int>
            {
                {'G', -15 }, {'P', -14 },{'L', -14 },{'K', -13 },{'Y', -5 },{'W', 13 },{'X',-1 }
            } },
            {'X', new Dictionary<char, int>
            {
                {'G', -1 }, {'P', -1 },{'L', -1 },{'K', -1 },{'Y', -1 },{'W',-1 },{'X',-1 }
            } },
        };

        //Smith–Waterman algorithm
        public static void AlignMatrix(MatrixAlignment MA, int gap, Dictionary<char, Dictionary<char, int>> scoreTable)
        {
            var _maxScore = 0;
            for (int j = 1; j < MA.QueryLength + 1; j++)
            {
                for (int i = 1; i < MA.ProteinLength + 1; i++)
                {
                    var lu = MA.matrix[j - 1][i-1].MaxCost + scoreTable[MA.Query[j-1]][MA.Protein[i-1]] > 0 ? MA.matrix[j - 1][i - 1].MaxCost + scoreTable[MA.Query[j - 1]][MA.Protein[i - 1]] : 0;
                    var l = MA.matrix[j][i -1].MaxCost - gap > 0 ? MA.matrix[j][i - 1].MaxCost - gap : 0;
                    var u = MA.matrix[j-1][i].MaxCost - gap > 0 ? MA.matrix[j - 1][i].MaxCost : 0;                   

                    // if same score, tend to choose lu > l > u.
                    if ( lu >= l && lu >= u)
                    {
                        MA.matrix[j][i].MaxCost = lu;
                        MA.matrix[j][i].Source = new Tuple<int, int, int>(j-1, i-1, 0);

                        if (lu >= _maxScore)
                        {
                            _maxScore = lu;
                            MA.MaxScore = lu;
                            MA.MaxScorePos = new Tuple<int, int>(j, i);
                        }
            
                    }
                    else if (l > lu && l >= u)
                    {
                        MA.matrix[j][i].MaxCost = l;
                        MA.matrix[j][i].Source = new Tuple<int, int, int>(j, i - 1, 1);
                    }
                    else if (u > lu && u > l)
                    {
                        MA.matrix[j][i].MaxCost = u;
                        MA.matrix[j][i].Source = new Tuple<int, int, int>(j - 1, i, -1);
                    }
                }
            }
        }

        //Tuple<int, int, int> (pos j, pos i, from 0:up-left, 1:left, -1:up)
        public static Stack<Tuple<int, int, int>> Traceback(MatrixAlignment MA)
        {
            Stack<Tuple<int, int, int>> path = new Stack<Tuple<int, int, int>>();

            if (MA.MaxScore > 0)
            {
                int j = MA.MaxScorePos.Item1;
                int i = MA.MaxScorePos.Item2;
                while (j > 0 && i > 0 && MA.matrix[j][i].MaxCost > 0)
                {
                    int preJ = MA.matrix[j][i].Source.Item1;
                    int preI = MA.matrix[j][i].Source.Item2;
                    int direction = MA.matrix[j][i].Source.Item3;
                    path.Push(new Tuple<int, int, int>(preJ, preI, direction));

                    j = preJ;
                    i = preI;
                }
            }

            return path;
        }

        //Needleman–Wunsch algorithm with retriction
        public static void AlignMatrix_UpLeft(MatrixAlignment MA, int gap, Dictionary<char, Dictionary<char, int>> scoreTable)
        {
            //initial score be minus infinite.
            var _maxScore = -1084;
            for (int j = 1; j < MA.QueryLength + 1; j++)
            {
                for (int i = 1; i < MA.ProteinLength + 1; i++)
                {
                    var lu = MA.matrix[j - 1][i - 1].MaxCost + scoreTable[MA.Query[j - 1]][MA.Protein[i - 1]];

                    MA.matrix[j][i].MaxCost = lu;
                    MA.matrix[j][i].Source = new Tuple<int, int, int>(j - 1, i - 1, 0);
                }
            }

            if (MA.ProteinLength >= MA.QueryLength)
            {
                for (int i = MA.QueryLength; i < MA.ProteinLength + 1; i++)
                {
                    if (MA.matrix[MA.QueryLength][i].MaxCost > _maxScore)
                    {
                        _maxScore = MA.matrix[MA.QueryLength][i].MaxCost;
                        MA.MaxScore = _maxScore;
                        MA.MaxScorePos = new Tuple<int, int>(MA.QueryLength, i);
                    }
                   
                }
            }
        }

        public static char[][] SequenceComparison(MatrixAlignment MA, Stack<Tuple<int, int, int>> path)
        {
            char[][] seqs = new char[2][];
            seqs[0] = new char[path.Count];
            seqs[1] = new char[path.Count];

            int ind = 0;
            foreach (var pos in path)
            {
                if (pos.Item3 == 0)
                {
                    seqs[0][ind] = MA.Query[pos.Item1];
                    seqs[1][ind] = MA.Protein[pos.Item2];
                }
                else if (pos.Item3 == 1)
                {
                    seqs[0][ind] = MA.Query[pos.Item1];
                    seqs[1][ind] = '-';
                }
                else if (pos.Item3 == -1)
                {
                    seqs[0][ind] = '-';
                    seqs[1][ind] = MA.Protein[pos.Item2];
                }
                ind++;
            }

            return seqs;
        }
    
        public static void PrintMatrix(MatrixAlignment MA, string filePath)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                string line = "\t\t";
                foreach (var item in MA.Protein)
                {
                    line += item + "\t";
                }
                output.WriteLine(line);

                int ind = 1;
                foreach (var array in MA.matrix)
                {
                    if (ind == 1)
                    {
                        line = "\t";
                    }
                    else if (ind >= 2)
                    {
                        line = MA.Query[ind-2] + "\t";
                    }
                

                    foreach (var a in array)
                    {
                        line += a.MaxCost + "\t";
                    }

                    output.WriteLine(line);

                    ind++;
                }
            }
        }
    }
}
