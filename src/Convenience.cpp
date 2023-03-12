/*
 *Some convenience functions
 */
#include <RcppCommon.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*
 *returns a distance matrix of Y
 */
//[[Rcpp::export]]
NumericMatrix dist(NumericMatrix y, double alpha) {
  int n = y.rows();
  int d = y.cols();
  int i, j, k;
  NumericMatrix mat(n, n);
  if ((d == 1) & (alpha == 1)) {
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        mat(i, j) = abs(y(i, 0) - y(j, 0));
        mat(j, i) = mat(i, j);
      }
    }
  }
  if ((d > 1) & (alpha == 1)) {
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        for (k = 0; k < d; k++) {
          mat(i, j) += (y(i, k) - y(j, k)) * (y(i, k) - y(j, k));
        }
        mat(i, j) = sqrt(mat(i, j));
        mat(j, i) = mat(i, j);
      }
    }
  }
  if ((d == 1) & (alpha != 1)) {
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        mat(i, j) = pow(abs(y(i, 0) - y(j, 0)), alpha);
        mat(j, i) = mat(i, j);
      }
    }
  }
  if ((d > 1) & (alpha != 1)) {
    double b = alpha / 2;
    for (i = 0; i < n; i++) {
      for (j = i; j < n; j++) {
        for (k = 0; k < d; k++) {
          mat(i, j) += (y(i, k) - y(j, k)) * (y(i, k) - y(j, k));
        }
        mat(i, j) = pow(mat(i, j), b);
        mat(j, i) = mat(i, j);
      }
    }
  }
  return (mat);
}

/*
 *returns a permutation which rearranges its first argument into ascending, it begins from 0
 */

IntegerVector ordered(NumericVector x)
{
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [& ](int i, int j)
  {
    return x[i] < x[j];
  });
  return idx;
}

/*
 *get the row max from a matrix
 */

NumericVector row_max(NumericMatrix m)
{
  int nrow = m.nrow();
  NumericVector max(nrow);
  for (int i = 0; i < nrow; i++)
  {
    max[i] = Rcpp::max(m(i, _));
  }
  return max;
}

/*
 *combination number, n choose r
 */

int nCr(int n, int r)
{
  if (r > n) return 0;
  if (r * 2 > n) r = n - r;
  if (r == 0) return 1;

  int result = n;
  for (int i = 2; i <= r; ++i)
  {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

/*
 *all combinations of subsets of size r choose from {0, 1, ..., n - 1}
 */

void combination(IntegerVector arr, int n, int r, int index, IntegerVector data, int i, IntegerMatrix &subset, int &k)
{
  if (index == r)
  {
    for (int j = 0; j < r; j++)
    {
      subset(k, j) = data[j];
    }
    k++;
    return;
  }
  if (i >= n)
  {
    return;
  }
  data[index] = arr[i];
  combination(arr, n, r, index + 1, data, i + 1, subset, k);
  combination(arr, n, r, index, data, i + 1, subset, k);
}

IntegerMatrix getallsubset(int n, int r)
{
  IntegerVector arr(n);
  for (int i = 0; i < n; i++)
  {
    arr[i] = i;
  }
  int L = nCr(n, r);
  IntegerMatrix subset(L, r);
  IntegerVector data(r);
  int k = 0;
  combination(arr, n, r, 0, data, 0, subset, k);
  return (subset);
}

/*
 *set difference on two vectors
 */

IntegerVector set_diff(IntegerVector set1, IntegerVector set2)
{
  int n = set1.size();
  int m = set2.size();
  IntegerVector index(n);
  IntegerVector diff(n - m);
  for (int i = 0; i < m; i++)
  {
    index[set2[i]] = 1;
  }
  int k = 0;
  for (int i = 0; i < n; i++)
  {
    if (index[i] == 0)
    {
      diff[k] = i;
      k++;
    }
  }
  return (diff);
}

/*
 *the same as which() in R, but the index begins from 0
 */

IntegerVector which(IntegerVector x, int y)
{
  int n = x.size();
  vector<int> a;
  for (int i = 0; i < n; i++)
  {
    if (x[i] == y)
    {
      a.push_back(i);
    }
  }
  return Rcpp::IntegerVector(a.begin(), a.end());
}

/*
 *initize the frame
 */

DataFrame init(int max_node)
{
  IntegerVector node(max_node);
  IntegerVector parent(max_node);
  IntegerVector leftChild(max_node);
  IntegerVector rightChild(max_node);
  NumericVector split(max_node);
  NumericVector statistic(max_node);
  NumericVector pval(max_node, NA_REAL);
  IntegerVector
    var (max_node);
  IntegerVector isleaf(max_node, -1);
  IntegerVector n(max_node);
  DataFrame frame = DataFrame::create(Named("node") = clone(node),
                                      Named("parent") = clone(parent),
                                      Named("leftChild") = clone(leftChild),
                                      Named("rightChild") = clone(rightChild),
                                      Named("statistic") = clone(statistic),
                                      Named("pval") = clone(pval),
                                      Named("split") = clone(split),
                                      Named("var") = clone(var),
                                      Named("isleaf") = clone(isleaf),
                                      Named("n") = clone(n));
  return (frame);
}
