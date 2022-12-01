\name{htt_control}
\alias{htt_control}

\title{
  Control for Hypothesis Testing Tree
}

\description{
  Various parameters that control aspects of the \code{\link{HTT}} function.
}

\usage{
htt_control(teststat = c("energy0", "energy1"),
             testtype = c("permutation", "fastpermutation"),
             alpha = 1, pt = 0.05, minsplit = 30,
             minbucket = round(minsplit/3),
             R = 199, nmin = 1000)
}

\arguments{
  \item{teststat}{ a character specifying the type of the test statistic to be applied.
                   It can be \code{teststat = "energy0"} or \code{teststat = "energy1"}.
                   Default is \code{teststat = "energy0"}.
}
  \item{testtype}{ a character specifying how to compute the distribution of the test statistic.
                   It can be \code{testtype = "permutation"} or \code{testtype = "fastpermutation"}.
                   For \code{testtype = "fastpermutation"}, it will not perform the permutation tests
                   on the node with more than \code{nmin} observations.
                   Default is \code{testtype = "permutation"}.
}
  \item{alpha}{ the exponent on Euclidean distance in (0,2] (for regression tree).
                Default is \code{alpha = 1}.
}
  \item{pt}{ the p-value of the permutation test must be less than in order to implement a split.
             If \code{pt = 1}, hypothesis testing tree will fully split
             without performing the permutation tests. Default is \code{pt = 0.05}.
}
  \item{minsplit}{ the minimum number of observations in a node
                   in order to be considered for splitting.
                   Default is \code{minsplit = 30}.
}
  \item{minbucket}{ the minimum number of observations in a terminal node.
                    Default is \code{minbucket = round(minsplit/3)}.
}
  \item{R}{ the number of permutation replications are used to simulated
            the distribution of the test statistic.
            Default is \code{R = 199}.
}
  \item{nmin}{ the minimum number of observations in a node that does not require
               the permutation test (for \code{testtype = "fastpermutation"}).
               Default is \code{nmin = 1000}.
}
}

\details{
  The arguments \code{teststat}, \code{testtype} and \code{pt} determine
  the hypothesis testing of each split.
  The argument \code{R} is the number of permutations to be used.
  For the dataset with more than 2000 observations, \code{testtype = "fastpermutation"} will be useful to save time.
}

\value{
  A list containing the options.
}

\seealso{
\code{\link{HTT}}, \code{\link{htt.object}}
}

\examples{
## choose the teststat as "energy1"
htt_control(teststat = "energy1")

## choose the p-value 0.01
htt_control(pt = 0.01)

## choose the alpha to 0.5
htt_control(alpha = 0.5)

## change the minimum number of observations in a terminal node
htt_control(minbucket = 7)

## reduce the number of permutation replications to save time
htt_control(R = 99)
}

\keyword{misc}
