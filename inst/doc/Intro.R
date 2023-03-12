## ---- echo=FALSE, message=FALSE, warning=FALSE--------------------------------
knitr::opts_chunk$set(comment = "#", warning = FALSE, eval = TRUE, message = FALSE)
set.seed(1)
library(HTT)

## ---- fig.width = 10, fig.height = 8, fig.align = "center", fig.retina = 2----
data("Boston", package = "MASS")
# set the p-value of the permutation test to 0.01
htt_boston <- HTT(medv ~ . , data = Boston, controls = htt_control(pt = 0.01))
htt_boston
# print the split information
htt_boston$frame
# Visualize HTT
plot(htt_boston)

## ---- fig.width = 9, fig.height = 7, fig.align = "center", fig.retina = 2-----
htt_iris <- HTT(Species ~., data = iris, controls = htt_control(pt = 0.01))
plot(htt_iris, layout = "tree")
# prediction 
table(predict(htt_iris), iris[, 5])

## -----------------------------------------------------------------------------
data("ENB")
set.seed(1)
idx = sample(1:nrow(ENB), floor(nrow(ENB)*0.8))
train = ENB[idx, ]
test = ENB[-idx, ]
htt_enb = HTT(cbind(Y1, Y2) ~ . , data = train, controls = htt_control(pt = 0.05, R = 99))
# prediction
pred = predict(htt_enb, newdata = test)
test_y = test[, 9:10]
# MAE
colMeans(abs(pred - test_y))
# MSE
colMeans(abs(pred - test_y)^2)

