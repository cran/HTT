\name{printsplit}
\alias{printsplit}

\title{
  Displays split table for Fitted htt Object
}

\description{
  Display the split table for fitted \code{htt} object.
}

\usage{
printsplit(object)
}

\arguments{
  \item{object}{ fitted model object of class \code{htt} returned
                 by the \code{HTT} function. }
}

\value{
  Display the split table.
}

\seealso{
\code{\link{HTT}}, \code{\link{htt.object}}
}

\examples{
irishtt = HTT(Species ~., data = iris)
printsplit(irishtt)
}
