using System;
using System.Collections.Generic;
using System.Text;

namespace Search_pool
{
    public class Target
    {
        public int Addup_Score { get; private set; }

        public string Protein { get; private set; }

        public string[] Matched_sequences { get; set; }

        public int[] SingleScores { get; set; }


        public static Target SetTarget(int score, string protein)
        {
            Target target = new Target();
            target.Addup_Score = score;
            target.Protein = protein;

            return target;
        }
    }
}
