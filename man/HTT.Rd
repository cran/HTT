\name{Hypothesis Testing Tree}
\alias{HTT}

\title{
Hypothesis Testing Tree
}

\description{
 Fit a hypothesis testing tree.
}

\usage{
HTT(X, Y, method, distance, controls = htt_control())
}

\arguments{
  \item{X}{ predictor variables in a vector, matrix or data.frame format. }
  \item{Y}{ response variables in a vector, factor or matrix format. }
  \item{method}{ \code{"regression"} or \code{"classification"}.
                 If \code{method} is missing then the routine tries to make an intelligent guess.
                 If \code{Y} is factor, then \code{method = "classification"}.
                 If \code{Y} is numeric vector or numeric matrix, then \code{method = "classification"}.
                 }
  \item{distance}{ If \code{distance} is missing, then Euclidean distance with exponent alpha is used for regression tree,
                   0-1 distance is used for classification tree.
                   Otherwise, use the \code{distance} as the distance matrix of \code{Y}. }
  \item{controls}{ a list of options that control details of the \code{HTT} algorithm. See \code{\link{htt_control}}. }
}

\details{
  Hypothesis testing trees examines the distribution difference over two child nodes by the binary partitioning in a hypothesis testing framework. At each split, it finds the maximum distribution difference over all possible binary partitions, the test statistic is based on generalized energy distance. The permutation test is used to estimate the p-value of the hypothesis testing.
}

\value{
  An object of class \code{htt}. See \code{\link{htt.object}}.
}

\author{
Jiaqi Hu
}

\examples{
## regression
data("Boston", package = "MASS")
Bostonhtt <- HTT(Boston[, 1:13], Boston[, 14], controls = htt_control(R = 99))
plot(Bostonhtt)
mean((Boston[, 14] - predict(Bostonhtt))^2)

## classification
irishtt <- HTT(iris[, 1:4], iris[, 5])
plot(irishtt)
mean(iris[, 5] == predict(irishtt))
}

\seealso{
  \code{\link{htt_control}}, \code{\link{print.htt}}, \code{\link{plot.htt}}, \code{\link{predict.htt}}
}

\keyword{HTT}
