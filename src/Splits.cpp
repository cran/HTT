/*
 *Binary splits for continuous/ordinal/nominal variables
 */
#include "Convenience.h"

/*Search for the optimal cutpoint in a ordered variable x maximizing the test statistic
 *\param x: as.numeric(x) (full sample)
 *\param sub_ind: the indices of the learning sample at the full sample
 *\param dmat: the distance matrix of response variable y (full sample)
 *\param var_type: 0: continuous, 1: ordinal
 *\param teststat: the type of the test statistic
 *\param minbucket: minimum number of sample in a leaf node
 *\param s: the sum of distance matrix of response variable y (learning sample)
 */

List split_ordered(NumericVector x, IntegerVector sub_ind, NumericMatrix dmat, int var_type,
                   string teststat, int minbucket, double s)
{
  double A = 0, B = 0, C = 0, maxstat = 0, temp = 0;
  int cutpoint = 0;
  int n = sub_ind.size();
  NumericVector sub_x(n);
  for (int i = 0; i < n; i++)
  {
    sub_x[i] = x(sub_ind[i]);
  }
  IntegerVector ord = ordered(sub_x);
  IntegerVector ind(n);
  for (int i = 0; i < n; i++)
  {
    ind[i] = sub_ind[ord[i]];
  }
  for (int i = 0; i < minbucket - 1; i++)
  {
    for (int j = 0; j < minbucket - 1; j++)
    {
      A += dmat(ind[j], ind[i]);
    }
  }
  for (int i = 0; i < minbucket - 1; i++)
  {
    for (int j = minbucket - 1; j < n; j++)
    {
      C += dmat(ind[j], ind[i]);
    }
  }
  B = s - A - 2 * C;
  if (teststat == "energy0")
  {
    for (int k = minbucket - 1; k < n - minbucket; k++)
    {
      for (int i = 0; i <= k; i++)
      {
        A += 2* dmat(ind[i], ind[k]);
      }
      for (int i = k; i < n; i++)
      {
        B -= 2* dmat(ind[i], ind[k]);
      }
      C = (s - A - B) / 2;
      temp = (2 *C / ((k + 1) *(n - k - 1)) - A / ((k + 1) *(k + 1)) - B / ((n - k - 1) *(n - k - 1))) *sqrt((k + 1) *(n - k - 1)) / log(log(n));
      if ((temp > maxstat) & (sub_x[ord[k]] != sub_x[ord[k + 1]]))
      {
        maxstat = temp;
        cutpoint = k;
      }
    }
  }
  else
  {
    for (int k = minbucket - 1; k < n - minbucket; k++)
    {
      for (int i = 0; i <= k; i++)
      {
        A += 2* dmat(ind[i], ind[k]);
      }
      for (int i = k; i < n; i++)
      {
        B -= 2* dmat(ind[i], ind[k]);
      }
      C = (s - A - B) / 2;
      temp = (2 *C / ((k + 1) *(n - k - 1)) - A / ((k + 1) *(k + 1)) - B / ((n - k - 1) *(n - k - 1))) *(k + 1) *(n - k - 1) / n;
      if ((temp > maxstat) & (sub_x[ord[k]] != sub_x[ord[k + 1]]))
      {
        maxstat = temp;
        cutpoint = k;
      }
    }
  }
  double split = (sub_x[ord[cutpoint]] + sub_x[ord[cutpoint + 1]]) / 2;
  if (var_type == 0)
  {
    // continuous
    return List::create(Named("maxstat") = maxstat, Named("cutpoint") = cutpoint + 1, Named("split") = split);
  }
  else
  {
    // ordinal
    NumericVector ux = sort_unique(sub_x);
    vector<int> left_leaf;
    vector<int> right_leaf;
    for (int i = 0; i < ux.size(); i++)
    {
      if (ux[i] < split)
      {
        left_leaf.push_back(ux[i]);
      }
      else
      {
        right_leaf.push_back(ux[i]);
      }
    }
    return List::create(Named("maxstat") = maxstat, Named("cutpoint") = cutpoint + 1, Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf);
  }
}

/*Search for the optimal cutpoint in an unordered factor x maximizing the test statistic
 *\param x: as.numeric(x) (full sample)
 *\param sub_ind: the indices of the learning sample at the full sample
 *\param dmat: the distance matrix of response variable y (full sample)
 *\param teststat: the type of the test statistic
 *\param minbucket: minimum number of sample in a leaf node
 *\param s: the sum of distance matrix of response variable y (learning sample)
 */

List split_unordered(NumericVector x, IntegerVector sub_ind, NumericMatrix dmat,
                     string teststat, int minbucket, double s)
{
  int n = sub_ind.size();
  IntegerVector sub_x(n);
  for (int i = 0; i < n; i++)
  {
    sub_x[i] = x(sub_ind[i]);
  }
  IntegerVector ux = sort_unique(sub_x);
  IntegerVector nn = table(sub_x);
  int m = nn.size(), cutpoint = 0;
  if (m == 1)
  {
    IntegerVector left_leaf;
    IntegerVector right_leaf = { 1 };
    return List::create(Named("maxstat") = 0, Named("cutpoint") = 0, Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf);
  }
  double temp = 0, maxstat = 0;
  IntegerVector left_index;
  NumericMatrix mat(m, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      IntegerVector indi = which(sub_x, ux[i]);
      IntegerVector indj = which(sub_x, ux[j]);
      int ni = indi.size();
      int nj = indj.size();
      for (int ii = 0; ii < ni; ii++)
      {
        for (int jj = 0; jj < nj; jj++)
        {
          mat(i, j) += dmat(sub_ind[indi[ii]], sub_ind[indj[jj]]);
        }
      }
      mat(j, i) = mat(i, j);
    }
  }
  IntegerVector arr(m);
  for (int i = 0; i < m; i++)
  {
    arr[i] = i;
  }
  if (teststat == "energy0")
  {
    for (int r = 1; r <= m / 2; r++)
    {
      int L = nCr(m, r);
      IntegerMatrix subset = getallsubset(m, r);
      for (int k = 0; k < L; k++)
      {
        IntegerVector ind1 = subset(k, _);
        IntegerVector ind2 = set_diff(arr, ind1);
        int n1 = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          n1 += nn[ind1[i]];
        }
        int n2 = n - n1;
        if ((n1 >= minbucket) & (n2 >= minbucket))
        {
          double A = 0, B = 0, C = 0;
          for (int i = 0; i < ind1.size(); i++)
          {
            for (int j = 0; j < ind1.size(); j++)
              A += mat(ind1[i], ind1[j]);
          }
          for (int i = 0; i < ind2.size(); i++)
          {
            for (int j = 0; j < ind2.size(); j++)
              B += mat(ind2[i], ind2[j]);
          }
          C = (s - A - B) / 2;
          temp = (2 *C / (n1 *n2) - A / (n1 *n1) - B / (n2 *n2)) *sqrt(n1 *n2) / log(log(n));
          if (temp > maxstat)
          {
            maxstat = temp;
            cutpoint = k;
            left_index = subset(cutpoint, _);
          }
        }
      }
    }
  }
  else
  {
    for (int r = 1; r <= m / 2; r++)
    {
      int L = nCr(m, r);
      IntegerMatrix subset = getallsubset(m, r);
      for (int k = 0; k < L; k++)
      {
        IntegerVector ind1 = subset(k, _);
        IntegerVector ind2 = set_diff(arr, ind1);
        int n1 = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          n1 += nn[ind1[i]];
        }
        int n2 = n - n1;
        if ((n1 >= minbucket) & (n2 >= minbucket))
        {
          double A = 0, B = 0, C = 0;
          for (int i = 0; i < ind1.size(); i++)
          {
            for (int j = 0; j < ind1.size(); j++)
              A += mat(ind1[i], ind1[j]);
          }
          for (int i = 0; i < ind2.size(); i++)
          {
            for (int j = 0; j < ind2.size(); j++)
              B += mat(ind2[i], ind2[j]);
          }
          C = (s - A - B) / 2;
          temp = (2 *C / (n1 *n2) - A / (n1 *n1) - B / (n2 *n2)) *n1 *n2 / n;
          if (temp > maxstat)
          {
            maxstat = temp;
            cutpoint = k;
            left_index = subset(cutpoint, _);
          }
        }
      }
    }
  }
  cutpoint = 0;
  for (int i = 0; i < left_index.size(); i++)
  {
    cutpoint += nn[left_index[i]];
  }
  int size_left = left_index.size();
  IntegerVector left_leaf(size_left);
  for (int i = 0; i < size_left; i++)
  {
    left_leaf[i] = ux[left_index[i]];
  }
  IntegerVector right_index = set_diff(arr, left_index);
  int size_right = right_index.size();
  IntegerVector right_leaf(size_right);
  for (int i = 0; i < size_right; i++)
  {
    right_leaf[i] = ux[right_index[i]];
  }
  if (min(left_leaf) < min(right_leaf))
  {
    return List::create(Named("maxstat") = maxstat, Named("cutpoint") = cutpoint, Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf);
  }
  else
  {
    return List::create(Named("maxstat") = maxstat, Named("cutpoint") = cutpoint, Named("left_leaf") = right_leaf, Named("right_leaf") = left_leaf);
  }
}

/*Greedy Search for the cutpoint in an unordered factor x maximizing the test statistic
 *\param x: as.numeric(x) (full sample)
 *\param sub_ind: the indices of the learning sample at the full sample
 *\param dmat: the distance matrix of response variable y (full sample)
 *\param teststat: the type of the test statistic
 *\param minbucket: minimum number of sample in a leaf node
 *\param s: the sum of distance matrix of response variable y (learning sample)
 */

List split_unordered_greedy(NumericVector x, IntegerVector sub_ind, NumericMatrix dmat,
                            string teststat, int minbucket, double s)
{
  int n = sub_ind.size();
  IntegerVector sub_x(n);
  for (int i = 0; i < n; i++)
  {
    sub_x[i] = x(sub_ind[i]);
  }
  IntegerVector ux = sort_unique(sub_x);
  IntegerVector nn = table(sub_x);
  int m = nn.size(), cutpoint = 0;
  if (m == 1)
  {
    IntegerVector left_leaf;
    IntegerVector right_leaf = { 1 };
    return List::create(Named("maxstat") = 0, Named("cutpoint") = 0, Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf);
  }
  double temp = 0, V = 0;
  NumericMatrix mat(m, m);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      IntegerVector indi = which(sub_x, ux[i]);
      IntegerVector indj = which(sub_x, ux[j]);
      int ni = indi.size();
      int nj = indj.size();
      for (int ii = 0; ii < ni; ii++)
      {
        for (int jj = 0; jj < nj; jj++)
        {
          mat(i, j) += dmat(sub_ind[indi[ii]], sub_ind[indj[jj]]);
        }
      }
      mat(j, i) = mat(i, j);
    }
  }
  IntegerVector arr(m);
  for (int i = 0; i < m; i++)
  {
    arr[i] = i;
  }
  vector<int> left_index;
  vector<int> right_index;
  for (int i = 0; i < m; i++)
  {
    right_index.push_back(i);
  }
  if (teststat == "energy0")
  {
    while (left_index.size() < m - 1)
    {
      double Vhat = 0;
      int bhat;
      for (int b: right_index)
      {
        left_index.push_back(b);
        IntegerVector ind1 = IntegerVector(left_index.begin(), left_index.end());
        IntegerVector ind2 = set_diff(arr, ind1);
        left_index.pop_back();
        int n1 = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          n1 += nn[ind1[i]];
        }
        int n2 = n - n1;
        double A = 0, B = 0, C = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          for (int j = 0; j < ind1.size(); j++)
            A += mat(ind1[i], ind1[j]);
        }
        for (int i = 0; i < ind2.size(); i++)
        {
          for (int j = 0; j < ind2.size(); j++)
            B += mat(ind2[i], ind2[j]);
        }
        C = (s - A - B) / 2;
        temp = (2 *C / (n1 *n2) - A / (n1 *n1) - B / (n2 *n2)) *sqrt(n1 *n2) / log(log(n));
        if (temp > Vhat)
        {
          Vhat = temp;
          bhat = b;
        }
      }
      int n1 = 0;
      for (int i = 0; i < left_index.size(); i++)
      {
        n1 += nn[left_index[i]];
      }
      int n2 = n - n1;
      if ((Vhat > V) | (n1 < minbucket) | (n2 < minbucket))
      {
        left_index.push_back(bhat);
        right_index.erase(std::find(right_index.begin(), right_index.end(), bhat));
        V = Vhat;
      }
      else
      {
        break;
      }
    }
  }
  else
  {
    while (left_index.size() < m - 1)
    {
      double Vhat = 0;
      int bhat;
      for (int b: right_index)
      {
        left_index.push_back(b);
        IntegerVector ind1 = IntegerVector(left_index.begin(), left_index.end());
        IntegerVector ind2 = set_diff(arr, ind1);
        left_index.pop_back();
        int n1 = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          n1 += nn[ind1[i]];
        }
        int n2 = n - n1;
        double A = 0, B = 0, C = 0;
        for (int i = 0; i < ind1.size(); i++)
        {
          for (int j = 0; j < ind1.size(); j++)
            A += mat(ind1[i], ind1[j]);
        }
        for (int i = 0; i < ind2.size(); i++)
        {
          for (int j = 0; j < ind2.size(); j++)
            B += mat(ind2[i], ind2[j]);
        }
        C = (s - A - B) / 2;
        temp = (2 *C / (n1 *n2) - A / (n1 *n1) - B / (n2 *n2)) *n1 *n2 / n;
        if (temp > Vhat)
        {
          Vhat = temp;
          bhat = b;
        }
      }
      int n1 = 0;
      for (int i = 0; i < left_index.size(); i++)
      {
        n1 += nn[left_index[i]];
      }
      int n2 = n - n1;
      if ((Vhat > V) | (n1 < minbucket) | (n2 < minbucket))
      {
        left_index.push_back(bhat);
        right_index.erase(std::find(right_index.begin(), right_index.end(), bhat));
        V = Vhat;
      }
      else
      {
        break;
      }
    }
  }
  for (int i = 0; i < left_index.size(); i++)
  {
    cutpoint += nn[left_index[i]];
  }
  int size_left = left_index.size();
  IntegerVector left_leaf(size_left);
  for (int i = 0; i < size_left; i++)
  {
    left_leaf[i] = ux[left_index[i]];
  }
  int size_right = right_index.size();
  IntegerVector right_leaf(size_right);
  for (int i = 0; i < size_right; i++)
  {
    right_leaf[i] = ux[right_index[i]];
  }
  left_leaf = sort_unique(left_leaf);
  right_leaf = sort_unique(right_leaf);
  if (min(left_leaf) < min(right_leaf))
  {
    return List::create(Named("maxstat") = V, Named("cutpoint") = cutpoint, Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf);
  }
  else
  {
    return List::create(Named("maxstat") = V, Named("cutpoint") = cutpoint, Named("left_leaf") = right_leaf, Named("right_leaf") = left_leaf);
  }
}

/*Search for the optimal cutpoint and variable by maximizing the test statistic
 *\param x: as.numeric(x) a numeric matrix
 *\param sub_ind: the indices of the learning sample at the full sample
 *\param dmat: the distance matrix of response variable y (full sample)
 *\param var_type: 0: continuous, 1: ordinal, 2: nominal(no more than 10 levels), 3: nominal(more than 10 levels)
 *\param teststat: the type of the test statistic
 *\param minbucket: minimum number of sample in a leaf node
 *\param s: the sum of distance matrix of response variable y (learning sample)
 */

List Split(int Node, NumericMatrix x, IntegerVector sub_ind, NumericMatrix dmat, IntegerVector var_type,
           string teststat, int minbucket, double s, DataFrame frame)
{
  int n = sub_ind.size(), p = x.cols();
  int optimal_dim = 0;
  double optimal_stat = 0, split = 0, maxstat = 0;
  IntegerVector left_leaf, right_leaf;
  for (int k = 0; k < p; k++)
  {
    if (var_type[k] < 2)
    {
      List res = split_ordered(x(_, k), sub_ind, dmat, var_type[k], teststat, minbucket, s);
      maxstat = res["maxstat"];
      if (maxstat > optimal_stat)
      {
        optimal_stat = maxstat;
        optimal_dim = k;
        if (var_type[k] == 0)
        {
          split = res["split"];
        }
        else
        {
          left_leaf = res["left_leaf"];
          right_leaf = res["right_leaf"];
        }
      }
    }
    else if (var_type[k] == 2)
    {
      NumericVector ux(sub_ind.size());
      for (int ii = 0; ii < sub_ind.size(); ii++)
      {
        ux[ii] = x(sub_ind[ii], k);
      }
      int ux_size = unique(ux).size();
      if (ux_size <= 10)
      {
        List res = split_unordered(x(_, k), sub_ind, dmat, teststat, minbucket, s);
        maxstat = res["maxstat"];
        if (maxstat > optimal_stat)
        {
          optimal_stat = maxstat;
          optimal_dim = k;
          left_leaf = res["left_leaf"];
          right_leaf = res["right_leaf"];
        }
      }
      else
      {
        List res = split_unordered_greedy(x(_, k), sub_ind, dmat, teststat, minbucket, s);
        maxstat = res["maxstat"];
        if (maxstat > optimal_stat)
        {
          optimal_stat = maxstat;
          optimal_dim = k;
          left_leaf = res["left_leaf"];
          right_leaf = res["right_leaf"];
        }
      }
    }
  }
  IntegerVector node = frame["node"];
  node[Node] = Node;
  IntegerVector nn = frame["n"];
  nn[Node] = n;
  NumericVector Split = frame["split"];
  NumericVector statistic = frame["statistic"];
  IntegerVector
    var = frame["var"];
  statistic[Node] = optimal_stat;
  if (optimal_stat > 0)
  {
    var[Node] = optimal_dim + 1;
    if (var_type[optimal_dim] == 0)
    {
      vector<int> left_index;
      vector<int> right_index;
      Split[Node] = split;
      for (int i = 0; i < n; i++)
      {
        if (x(sub_ind[i], optimal_dim) < split)
        {
          left_index.push_back(sub_ind[i]);
        }
        else
        {
          right_index.push_back(sub_ind[i]);
        }
      }
      return List::create(Named("optimal_stat") = optimal_stat, Named("optimal_dim") = optimal_dim,
                          Named("split") = split,
                          Named("left_index") = left_index, Named("right_index") = right_index);
    }
    else
    {
      vector<int> left_index;
      vector<int> right_index;
      for (int i = 0; i < n; i++)
      {
        if (which(left_leaf, x(sub_ind[i], optimal_dim)).size() > 0)
        {
          left_index.push_back(sub_ind[i]);
        }
        else
        {
          right_index.push_back(sub_ind[i]);
        }
      }
      Split[Node] = NA_REAL;
      return List::create(Named("optimal_stat") = optimal_stat, Named("optimal_dim") = optimal_dim,
                          Named("left_leaf") = left_leaf, Named("right_leaf") = right_leaf,
                          Named("left_index") = left_index, Named("right_index") = right_index);
    }
  }
  else
  {
    var[Node] = NA_INTEGER;
    IntegerVector isleaf = frame["isleaf"];
    isleaf[Node] = 1;
    return List::create(Named("optimal_stat") = optimal_stat, Named("optimal_dim") = NA_INTEGER, Named("split") = NA_REAL);
  }
}

/*
 *update the frame information
 */

void update(NumericMatrix X, NumericMatrix dmat, NumericMatrix perm_stat_mat, DataFrame frame, IntegerVector var_type,
            string teststat, int minbucket, IntegerVector &where, int &flag, double pt)
{
  IntegerVector node = frame["node"];
  IntegerVector parent = frame["parent"];
  IntegerVector leftChild = frame["leftChild"];
  IntegerVector rightChild = frame["rightChild"];
  NumericVector split = frame["split"];
  NumericVector statistic = frame["statistic"];
  NumericVector pval = frame["pval"];
  IntegerVector isleaf = frame["isleaf"];
  IntegerVector var = frame["var"];
  IntegerVector n = frame["n"];
  int R = perm_stat_mat.rows();
  int current_node = which_max(node);
  IntegerVector current_leaf = which(isleaf[Rcpp::Range(0, current_node)], -1);
  if (current_leaf.size() == 0)
  {
    flag = 0;
    return;
  }
  NumericVector stat_leaf = statistic[current_leaf];
  int max_leaf = which_max(stat_leaf);
  max_leaf = current_leaf[max_leaf];
  NumericMatrix perm_stat(R, current_leaf.size());
  for (int i = 0; i < current_leaf.size(); i++)
  {
    perm_stat(_, i) = perm_stat_mat(_, current_leaf[i]);
  }
  NumericVector perm_stat1 = row_max(perm_stat);
  double pvalue = mean(perm_stat1 >= statistic[max_leaf]) / (R + 1) *R;
  pval[max_leaf] = pvalue;
  if (pvalue <= pt)
  {
    flag = 1;
    isleaf[max_leaf] = 0;
    parent[current_node + 1] = max_leaf;
    parent[current_node + 2] = max_leaf;
    leftChild[max_leaf] = current_node + 1;
    rightChild[max_leaf] = current_node + 2;
    IntegerVector sub_ind = which(where, max_leaf);
    int N = sub_ind.size();
    double s = 0;
    for (int i = 1; i < N; i++)
    {
      for (int j = 0; j < i; j++)
      {
        s += dmat(sub_ind[i], sub_ind[j]);
      }
    }
    s = s * 2;
    List res = Split(max_leaf, X, sub_ind, dmat, var_type, teststat, minbucket, s, frame);
    vector<int> left_index = res["left_index"];
    vector<int> right_index = res["right_index"];
    int nL = left_index.size();
    int nR = right_index.size();

    for (int i = 0; i < nL; i++)
    {
      where[left_index[i]] = current_node + 1;
    }
    for (int i = 0; i < nR; i++)
    {
      where[right_index[i]] = current_node + 2;
    }
  }
  else
  {
    flag = 0;
    isleaf[current_leaf] = 1;
  }
}

/*
 *update the frame information when no split
 */

void update1(int Node, int N, DataFrame frame)
{
  IntegerVector node = frame["node"];
  IntegerVector var = frame["var"];
  IntegerVector n = frame["n"];
  node[Node] = Node;
  var[Node] = NA_INTEGER;
  n[Node] = N;
}
