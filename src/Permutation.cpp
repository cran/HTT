/*
 *Binary splits
 */
#include "Convenience.h"
#include "Splits.h"

/*Search for the optimal cutpoint and variable by maximizing a two-sample statistic
 *\param sub_x: as.numeric(x)
 *\param sub_ind: the indices of the learning sample at the full sample
 *\param dmat: the distance matrix of response variable y (full sample)
 *\param var_type: 0: continuous, 1: ordinal, 2: nominal
 *\param teststat: the type of the test statistic
 *\param minbucket: minimum number of sample in a leaf node
 *\param s: the sum of distance matrix of response variable y (learning sample)
 */

Function set_seed("set.seed");

void Permute(int node, NumericMatrix x, IntegerVector sub_ind, NumericMatrix dmat,
             IntegerVector var_type, string teststat, int minbucket,
             double s, int R, NumericMatrix &perm_stat_mat)
{
  int n = sub_ind.size(), p = x.cols();
  double maxstat = 0, temp = 0;
  NumericVector perm_stat(R);
  for (int r = 0; r < R; r++)
  {
    maxstat = 0;
    set_seed(R *node + r);
    IntegerVector ind_perm = sample(sub_ind, sub_ind.size(), false);
    NumericMatrix xx = clone(x);
    for (int i = 1; i < n; i++)
    {
      xx(sub_ind[i], _) = x(ind_perm[i], _);
    }
    for (int k = 0; k < p; k++)
    {
      if (var_type[k] < 2)
      {
        List res = split_ordered(xx(_, k), sub_ind, dmat, var_type[k], teststat, minbucket, s);
        temp = res["maxstat"];
        if (temp > maxstat)
          maxstat = temp;
      }
      else if (var_type[k] == 2)
      {
        NumericVector ux(sub_ind.size());
        for (int ii = 0; ii < sub_ind.size(); ii++)
        {
          ux[ii] = xx(sub_ind[ii], k);
        }
        int ux_size = unique(ux).size();
        if (ux_size <= 10)
        {
          List res = split_unordered(xx(_, k), sub_ind, dmat, teststat, minbucket, s);
          temp = res["maxstat"];
          if (temp > maxstat)
            maxstat = temp;
        }
        else
        {
          List res = split_unordered_greedy(xx(_, k), sub_ind, dmat, teststat, minbucket, s);
          temp = res["maxstat"];
          if (temp > maxstat)
            maxstat = temp;
        }
      }
    }
    perm_stat[r] = maxstat;
  }
  perm_stat_mat(_, node) = perm_stat;
}
