#include <RcppCommon.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

NumericMatrix distance(NumericMatrix y, double alpha);
IntegerVector ordered(NumericVector x);
NumericVector row_max(NumericMatrix m);
int nCr(int n, int r);
void combination(IntegerVector arr, int n, int r, int index, IntegerVector data, int i, IntegerMatrix &subset, int &k);
IntegerMatrix getallsubset(int n, int r);
IntegerVector set_diff(IntegerVector set1, IntegerVector set2);
IntegerVector which(IntegerVector x, int y);
DataFrame init(int max_node);
