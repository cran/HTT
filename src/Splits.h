List split_ordered(NumericVector sub_x, IntegerVector sub_ind, NumericMatrix dmat, int var_type,
                   std::string teststat, int minbucket, double s);
List split_unordered(NumericVector sub_x, IntegerVector sub_ind, NumericMatrix dmat,
                     std::string teststat, int minbucket, double s);
List Split(int Node, NumericMatrix x, IntegerVector sub_ind, NumericMatrix dmat, IntegerVector var_type,
           std::string teststat, int minbucket, double s, DataFrame frame);
List split_unordered_greedy(NumericVector x, IntegerVector sub_ind, NumericMatrix dmat,
                            std::string teststat, int minbucket, double s);
void update(NumericMatrix X, NumericMatrix dmat, NumericMatrix perm_stat_mat, DataFrame frame, IntegerVector var_type,
            std::string teststat, int minbucket, IntegerVector &where, int &flag, double pt);
void update1(int Node, int N, DataFrame frame);
