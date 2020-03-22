using System;
using System.Collections.Generic;
using System.Text;
using Meta.Numerics.Functions;

namespace Search_pool
{
    public class EValue
    {
        /* 
        1.) gap existence penalty (INT2_MAX denotes infinite).
        2.) gap extension penalty (INT2_MAX denotes infinite).
        3.) decline to align penalty (this field is ignored)
        4.) Lambda
        5.) K
        6.) H
        7.) alpha
        8.) beta
        9.) C
        10.) alpha_v
        11.) sigma
        */
        public static double[][] pam30_values = new double[11][] {
            new double[11]{32767, 32767, 32767, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818 },
            new double[11]{7, 2, 32767, 0.305, 0.15, 0.87, 0.35, -3, 0.479087, 1.014010, 1.162730},
            new double[11]{6, 2, 32767, 0.287, 0.11, 0.68, 0.42, -4, 0.499980, 1.688060, 1.951430},
            new double[11]{5, 2, 32767, 0.264, 0.079, 0.45, 0.59, -7, 0.533009, 3.377010, 3.871950},
            new double[11]{10, 1, 32767, 0.309, 0.15, 0.88, 0.35, -3, 0.474741, 1.372050, 1.788770},
            new double[11]{9, 1, 32767, 0.294, 0.11, 0.61, 0.48, -6, 0.492716, 2.463920, 3.186150},
            new double[11]{8, 1, 32767, 0.270, 0.072, 0.40, 0.68, -10, 0.521286, 5.368130, 6.763480},
            new double[11]{15, 3, 32767, 0.339, 0.28, 1.70, 0.20, -0.5, 0.437688, 0.157089, 0.155299},
            new double[11]{14, 2, 32767, 0.337, 0.27, 1.62, 0.21, -0.8, 0.440010, 0.206970, 0.198524},
            new double[11]{14, 1, 32767, 0.333, 0.27, 1.43, 0.23, -1.4, 0.444817, 0.436301, 0.361947},
            new double[11]{13, 3, 32767, 0.338, 0.27, 1.69, 0.20, -0.5, 0.439086, 0.178973, 0.175436},
        };

        public static Blast_KarlinBlk Set_Blast_KarlinBlk(int ind, double paramC = 0)
        {
            return new Blast_KarlinBlk(pam30_values[ind][3], pam30_values[ind][4], Math.Log(pam30_values[ind][4]), pam30_values[ind][5], paramC);
        }

        public static Blast_GumbelBlk Set_Blast_GumbelBlk(int gapOpen, int gapExtent, int ind, int? db_length = null)
        {
            var g = (double)(gapOpen + gapExtent);
            var b = 2.0 * g * (pam30_values[0][6] - pam30_values[ind][6]);
            var beta = 2.0 * g * (pam30_values[0][9] - pam30_values[ind][9]);
            var tau = 2.0 * g * (pam30_values[0][9] - pam30_values[ind][10]);
            return new Blast_GumbelBlk(pam30_values[ind][3], pam30_values[ind][8], g, pam30_values[ind][6], 
                pam30_values[ind][9], pam30_values[ind][10], pam30_values[0][6], pam30_values[0][9], b, beta, tau, db_length, true);
        }

        public static double BLAST_KarlinStoE_simple(int score, Blast_KarlinBlk kbp, int searchsp)
        {
            var lambda = kbp.Lambda;
            var k = kbp.K;
            var h = kbp.H;

            if (lambda < 0.0 || k < 0.0 || h < 0.0)
            {
                return -1;
            }

            return searchsp * Math.Exp(-lambda * score + kbp.logK );
        }

        public static double BLAST_SpougeStoE(int y_, Blast_KarlinBlk kbp, Blast_GumbelBlk gbp, int m_, int n_)
        {
            // the score and lambda may have been rescaled.  We will compute the scaling factor
            // and use it to scale a, alpha, and Sigma similarly.
            double scale_factor = kbp.Lambda / gbp.Lambda;

            // the pair-wise e-value must be scaled back to db-wise e-value
            double db_scale_factor = (gbp.db_length.HasValue) ?
                    (double)gbp.db_length.Value / (double)n_ : 1.0;

            double lambda_ = kbp.Lambda;
            double k_ = kbp.K;
            double ai_hat_ = gbp.a * scale_factor;
            double bi_hat_ = gbp.b;
            double alphai_hat_ = gbp.Alpha * scale_factor;
            double betai_hat_ = gbp.Beta;
            double sigma_hat_ = gbp.Sigma * scale_factor;
            double tau_hat_ = gbp.Tau;

            /* here we consider symmetric matrix only */
            double aj_hat_ = ai_hat_;
            double bj_hat_ = bi_hat_;
            double alphaj_hat_ = alphai_hat_;
            double betaj_hat_ = betai_hat_;

            /* this is 1/sqrt(2.0*PI) */
            double const_val = 0.39894228040143267793994605993438;

            double m_li_y, vi_y, sqrt_vi_y, m_F, P_m_F;
            double n_lj_y, vj_y, sqrt_vj_y, n_F, P_n_F;
            double c_y, p1, p2, area;
            double e_value;

            m_li_y = m_ - (ai_hat_ * y_ + bi_hat_);
            vi_y = Math.Max(2.0 * alphai_hat_ / lambda_, alphai_hat_ * y_ + betai_hat_);
            sqrt_vi_y = Math.Sqrt(vi_y);
            m_F = m_li_y / sqrt_vi_y;
            P_m_F = AdvancedMath.Erfc(-m_F / Math.Sqrt(2.0)) / 2.0;
            p1 = m_li_y * P_m_F + sqrt_vi_y * const_val * Math.Exp(-0.5 * m_F * m_F);

            n_lj_y = n_ - (aj_hat_ * y_ + bj_hat_);
            vj_y = Math.Max(2.0 * alphaj_hat_ / lambda_, alphaj_hat_ * y_ + betaj_hat_);
            sqrt_vj_y = Math.Sqrt(vj_y);
            n_F = n_lj_y / sqrt_vj_y;
            P_n_F = AdvancedMath.Erfc(-n_F / Math.Sqrt(2.0)) / 2.0;
            p2 = n_lj_y * P_n_F + sqrt_vj_y * const_val * Math.Exp(-0.5 * n_F * n_F);

            c_y = Math.Max(2.0 * sigma_hat_ / lambda_, sigma_hat_ * y_ + tau_hat_);
            area = p1 * p2 + c_y * P_m_F * P_n_F;

            e_value = area * k_ * Math.Exp(-lambda_ * y_) * db_scale_factor;
            //ASSERT(e_value >= 0.0);

            return e_value;
        }

    }

    public struct Blast_KarlinBlk
    {
        public Blast_KarlinBlk(double lambda, double k, double lk, double h, double pc)
        {
            Lambda = lambda;
            K = k;
            logK = lk;
            H = h;
            paramC = pc;
        }
        public double Lambda { get; } /**< Lambda value used in statistics */
        public double K { get; } /**< K value used in statistics */
        public double logK { get; } /**< natural log of K value used in statistics */
        public double H { get; } /**< H value used in statistics */
        public double paramC { get; }  /**< for use in seed. */
    }

    public struct Blast_GumbelBlk
    {
        public Blast_GumbelBlk(double lambda, double c, double g, double _a, double alpha, double sigma, double _a_un, double alpha_un, double _b, double beta, double tau, int? _db, bool f)
        {
            Lambda = lambda;
            C = c;
            G = g;
            a = _a;
            Alpha = alpha;
            Sigma = sigma;
            a_un = _a_un;
            Alpha_un = alpha_un;
            b = _b;
            Beta = beta;
            Tau = tau;
            db_length = _db;
            filled = f;
        }

        public double Lambda { get; }   /**< the unscaled Lambda value */
        public double C { get; }
        public double G { get; }        /**< G is the total penalty for extension */
        public double a { get; }     /**< avg(L) = a     y + b    */
        public double Alpha { get; }     /**< var(L) = alpha y + beta */
        public double Sigma { get; }     /**< cov(L) = sigma y + tau  */
        public double a_un { get; }      /**< Ungapped a */
        public double Alpha_un { get; }  /**< Ungapped alpha */

        public double b { get; }         /**< 2*G*(a_un - a) */
        public double Beta { get; }      /**< 2*G*(alpha_un - alpha) */
        public double Tau { get; }       /**< 2*G*(alpha_un - Sigma) */

        public int? db_length { get; }    /**< total length of database */

        public bool filled { get; }    /**< flag indicate the values of gbp are prepared */

    }
}
