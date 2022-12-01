\name{plot.htt}
\alias{plot.htt}

\title{
  Plot an htt Object
}

\description{
  Visualize a \code{htt} object, several arguments can be passed to control the color and shape.
}

\usage{
\method{plot}{htt}(x, digits = 3,
    line.color = "blue",
    node.color = "black",
    line.type = c("straight", "curved"),
    layout = c("tree", "dendrogram"), \dots)
}

\arguments{
  \item{x}{ fitted model object of class \code{htt} returned
            by the \code{HTT} function.}
  \item{digits}{ the number of significant digits in displayed numbers.
                 Default is \code{digits = 3}.}
  \item{line.color}{ a character specifying the edge color.
                     Default is \code{line.color = "blue"}.}
  \item{node.color}{ a character specifying the node color.
                     Default is \code{node.color = "black"}.}
  \item{line.type}{ a character specifying the type of edge,
                    \code{line.type = "straight"} or \code{line.type = "curved"}.
                    Default is \code{line.type = "straight"}.}
  \item{layout}{ a character specifying the layout,
                 \code{layout = "tree"} or \code{layout = "dendrogram"}.
                 Default is \code{layout = "tree"}.}
  \item{\dots}{ additional print arguments.}
}

\details{
  This function is a method for the generic function \code{\link{plot}}, for objects of class \code{htt}.
}

\value{
  Visualize the hypothesis testing tree.
}

\seealso{
  \code{\link{print.htt}}, \code{\link{printsplit}}, \code{\link{predict.htt}}
}

\examples{
fit = HTT(iris[, 1:4], iris[, 5])
plot(fit)

# change the line color and node color
plot(fit, line.color = "black", node.color = "blue")

# change the line type
plot(fit, line.type = "curved")

# change the layout
plot(fit, layout = "dendrogram")
}
