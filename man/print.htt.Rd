\name{print.htt}
\alias{print.htt}

\title{
  Print a Fitted htt Object
}

\description{
  This function prints a \code{\link{htt.object}}.
  It is a method for the generic function \code{\link{print}} of class \code{htt}.
  It can be invoked by calling \code{print} for an object of the appropriate class,
  or directly by calling \code{print.htt} regardless of the class of the object.
}

\usage{
\method{print}{htt}(x, ...)
}

\arguments{
  \item{x}{ fitted model object of class \code{htt} returned
            by the \code{HTT} function. }
  \item{\dots}{ additional print arguments.}
}

\details{
  A semi-graphical layout of the contents of x$frame is printed.
  Indentation is used to convey the tree topology.
  Information for each node includes
  the node number,split rule, size and p-value.
  For the "class" method, the class probabilities are also printed.
}

\value{
  Visualize the hypothesis testing tree in a semi-graphical layout.
}

\seealso{
  \code{\link{htt.object}}, \code{\link{printsplit}}
}

\examples{
irishtt = HTT(Species ~., data = iris)
print(irishtt)
}
