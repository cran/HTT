\name{predict.htt}
\alias{predict.htt}

\title{
Predictions from a Fitted htt Object
}

\description{
Compute predictions from \code{htt} object.
}

\usage{
\method{predict}{htt}(object, newdata,
        type = c("response", "prob", "node"),
        ...)
}

\arguments{
  \item{object}{ fitted model object of class \code{htt} returned
                 by the \code{HTT} function. }
  \item{newdata}{ an optional data frame in which to look for variables
                  with which to predict, if omitted, the fitted values are used. }
  \item{type}{ a character string denoting the type of predicted value returned.
               For \code{type = "response"}, the mean of a numeric response and
               the predicted class for a categorical response is returned.
               For \code{type = "prob"} the matrix of class probabilities
               is returned for a categorical response. \code{type = "node"}
               returns an integer vector of terminal node identifiers. }
  \item{\dots}{ additional print arguments.}
}

\details{
  This function is a method for the generic function \code{predict}
  for class \code{htt}. It can be invoked by calling \code{predict} for an
  \code{object} of the appropriate class, or directly by calling \code{predict.htt}
  regardless of the class of the object.
}

\value{
  A list of predictions, possibly simplified to a numeric vector,
  numeric matrix or factor.

  If \code{type = "response"}:\cr
  the mean of a numeric response and
  the predicted class for a categorical response is returned.

  If \code{type = "prob"}:\cr
  the matrix of class probabilities
  is returned for a categorical response.

  If \code{type = "node"}:\cr
  an integer vector of terminal node identifiers is returned.
}

\seealso{
\code{\link{predict}}, \code{\link{htt.object}}
}

\examples{
irishtt <- HTT(Species ~., data = iris)

## the predicted class
predict(irishtt, type = "response")

## class probabilities
predict(irishtt, type = "prob")

## terminal node identifiers
predict(irishtt, type = "node")
}
