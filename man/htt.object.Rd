\name{htt.object}
\alias{htt.object}

\title{Hypothesis Testing Tree Object}

\description{
  A class for representing hypothesis testing tree.
}

\value{
\item{frame}{
  a dataframe about the split information. It contains following information:

  \code{node}: the node numbers in a split order.

  \code{parent}: the parent node number.

  \code{leftChild}: the left daughter node number, \code{NA} represents leaf.

  \code{rightChild}: the right daughter node number, \code{NA} represents leaf.

  \code{statistic}: the maximum test statistic of all possible splits within the node.

  \code{split}: the rule of the split. It is numeric for continuous covariate split,
  it is a character for non-numeric covariate split, the levels of two child nodes are
  stored in two braces.

  \code{pval}: approximate p-values estimated from permutation test.

  \code{isleaf}: 1 denotes terminal node and 0 denotes internal node.

  \code{n}: the number of observations reaching the node.

  \code{var}: the names of the variables used in the
  split at each node (leaf nodes are denoted by the label \code{"<leaf>"}).

  \code{yval}: the fitted value of the response at the node, if the dimension of response
  is larger than 1, it will presents as yval1, yval2, ... .

  \code{prob}: the probability of each class at the node, only visible for classification tree.
}

\item{where}{
  an integer vector of the same length as the number of observations in the
  root node, containing the row number of \code{frame} corresponding to
  the leaf node that each observation falls into.
}

\item{method}{
  the method used to grow the hypothesis testing tree, \code{"regression"} or \code{"classification"}.
}

\item{control}{
  a list of options that control the \code{HTT} algorithm. See \code{\link{htt_control}}.
}

\item{X}{
  a copy of the input \code{X} in a dataframe format.
}

\item{var.type}{
  a vector recording for each variables, 0 represents continuous,
  1 represents ordinal and 2 represents nominal variables.
}
}

\seealso{
  \code{\link{HTT}}, \code{\link{plot.htt}}, \code{\link{print.htt}}, \code{\link{predict.htt}}
}

\keyword{classes}
