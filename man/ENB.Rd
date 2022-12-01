\name{ENB}
\alias{ENB}
\docType{data}
\title{
Energy efficiency dataset
}
\description{
The data is about energy performance of buildings, containing eight input variables: relative compactness, surface area, wall area, roof area, overall height, orientation, glazing area, glazing area distribution and two output variables: heating load (HL) and cooling load (CL) of residential buildings. The goal is to predict two real valued responses from eight input variables. It can also be used as a multi-class classification problem if the response is rounded to the nearest integer.
}
\usage{data("ENB")}
\format{
  A data frame with 768 observations on the following 10 variables.
  \describe{
    \item{\code{X1}}{Relative Compactness}
    \item{\code{X2}}{Surface Area}
    \item{\code{X3}}{Wall Area}
    \item{\code{X4}}{Roof Area}
    \item{\code{X5}}{Overall Height}
    \item{\code{X6}}{Orientation}
    \item{\code{X7}}{Glazing Area}
    \item{\code{X8}}{Glazing Area Distribution}
    \item{\code{Y1}}{Heating Load}
    \item{\code{Y2}}{Cooling Load}
  }
}
\source{
UCI Machine Learning Repository: \url{https://archive.ics.uci.edu/ml/datasets/Energy+efficiency}.
}
\references{
A. Tsanas, A. Xifara: 'Accurate quantitative estimation of energy performance of residential buildings using statistical machine learning tools', Energy and Buildings, Vol. 49, pp. 560-567, 2012
}
\examples{
data(ENB)
set.seed(1)
train = sample(1:nrow(ENB), floor(nrow(ENB)*0.8))
train_x = ENB[train, 1:8]
train_y = ENB[train, 9:10]
test_x = ENB[-train, 1:8]
test_y = ENB[-train, 9:10]
htt_enb = HTT(train_x, train_y, controls = htt_control(pt = 0.05, R = 99))
# prediction
pred = predict(htt_enb, newdata = test_x)
# MAE
colMeans(abs(pred - test_y))
# MSE
colMeans(abs(pred - test_y)^2)
}
\keyword{datasets}
