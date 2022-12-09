/*
 *Growing a hypothesis test tree
 */
#include "Convenience.h"
#include "Splits.h"
#include "Permutation.h"

//[[Rcpp::export]]

List TreeGrow(NumericMatrix X, NumericMatrix dmat, IntegerVector var_type, List controls = R_NilValue)
{
  std::string teststat = controls["teststat"];
  std::string testtype = controls["testtype"];
  int R = controls["R"];
  int minsplit = controls["minsplit"];
  int minbucket = controls["minbucket"];
  int max_node = controls["maxnode"];
  int nmin = controls["nmin"];
  double pt = controls["pt"];
  DataFrame frame = init(max_node);
  IntegerVector where(X.rows());
  NumericMatrix perm_stat_mat(R, max_node);
  int flag = 1;
  int node = 0;
  IntegerVector sub_ind = which(where, node);
  int N = sub_ind.size();
  double s = sum(dmat);
  Split(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, frame);
  if ((testtype == "fastpermutation") & (N >= nmin)){
    for(int r = 0; r < R; r++){
      perm_stat_mat(r, node) = 0;
    }
  } else {
    Permute(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, R, perm_stat_mat);
  }
  while (flag & (node < max_node - 2))
  {
    update(X, dmat, perm_stat_mat, frame, var_type, teststat, minbucket, where, flag, pt);
    if (flag == 1)
    {
      node = node + 1;
      IntegerVector sub_ind = which(where, node);
      N = sub_ind.size();
      if (N >= minsplit)
      {
        double s = 0;
        for (int i = 1; i < N; i++)
        {
          for (int j = 0; j < i; j++)
          {
            s += dmat(sub_ind[i], sub_ind[j]);
          }
        }
        s = s * 2;
        Split(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, frame);
        if ((testtype == "fastpermutation") & (N >= nmin)){
          for(int r = 0; r < R; r++){
            perm_stat_mat(r, node) = 0;
          }
        } else {
          Permute(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, R, perm_stat_mat);
        }
      }
      else
      {
        update1(node, N, frame);
      }
      node = node + 1;
      sub_ind = which(where, node);
      N = sub_ind.size();
      if (N >= minsplit)
      {
        double s = 0;
        for (int i = 1; i < N; i++)
        {
          for (int j = 0; j < i; j++)
          {
            s += dmat(sub_ind[i], sub_ind[j]);
          }
        }
        s = s * 2;
        Split(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, frame);
        if ((testtype == "fastpermutation") & (N >= nmin)){
          for(int r = 0; r < R; r++){
            perm_stat_mat(r, node) = 0;
          }
        } else {
          Permute(node, X, sub_ind, dmat, var_type, teststat, minbucket, s, R, perm_stat_mat);
        }
      }
      else
      {
        update1(node, N, frame);
      }
    }
  }
  return List::create(Named("frame") = frame, Named("where") = where, Named("perm_stat_mat") = perm_stat_mat);
}
