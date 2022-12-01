# Hypothesis Testing Tree

## Regression Tree

    data("Boston", package = "MASS")
    # set the p-value of the permutation test to 0.01
    htt_boston <- HTT(X = Boston[, 1:13], Y = Boston[, 14], controls = htt_control(pt = 0.01))
    htt_boston

    #      Hypothesis Testing Tree 
    # 
    # node, split, n, pvalue
    # * denotes terminal node
    # 
    # [1] root   (n = 506, pvalue = 0)
    # |  [2] rm<=7.437   (n = 476, pvalue = 0)
    # |  |  [4] lstat<=15   (n = 314, pvalue = 0)
    # |  |  |  [6] rm<=6.797   (n = 256, pvalue = 0)
    # |  |  |  |  [8] lstat<=4.615   (n = 10) *
    # |  |  |  |  [9] lstat>4.615   (n = 246, pvalue = 0)
    # |  |  |  |  |  [12] rm<=6.543   (n = 212, pvalue = 0)
    # |  |  |  |  |  |  [14] lstat<=7.57   (n = 42) *
    # |  |  |  |  |  |  [15] lstat>7.57   (n = 170) *
    # |  |  |  |  |  [13] rm>6.543   (n = 34) *
    # |  |  |  [7] rm>6.797   (n = 58) *
    # |  |  [5] lstat>15   (n = 162, pvalue = 0)
    # |  |  |  [10] crim<=0.65402   (n = 46) *
    # |  |  |  [11] crim>0.65402   (n = 116, pvalue = 0)
    # |  |  |  |  [16] crim<=11.36915   (n = 77) *
    # |  |  |  |  [17] crim>11.36915   (n = 39) *
    # |  [3] rm>7.437   (n = 30) *

    # print the split information
    htt_boston$frame

    #    node parent leftChild rightChild  statistic pval    split     var isleaf   n
    # 1     1      0         2          3 2258.92680 0.00    7.437      rm      0 506
    # 2     2      1         4          5 1126.14057 0.00       15   lstat      0 476
    # 3     3      1        NA         NA   54.73540   NA   <leaf> ptratio      1  30
    # 4     4      2         6          7  750.08329 0.00    6.797      rm      0 314
    # 5     5      2        10         11  201.23810 0.00  0.65402    crim      0 162
    # 6     6      4         8          9  284.52923 0.00    4.615   lstat      0 256
    # 7     7      4        NA         NA   54.33706   NA   <leaf>   lstat      1  58
    # 8     8      6        NA         NA    0.00000   NA   <leaf>    <NA>      1  10
    # 9     9      6        12         13  188.93990 0.00    6.543      rm      0 246
    # 10   10      5        NA         NA   73.70296   NA   <leaf>     dis      1  46
    # 11   11      5        16         17  115.47482 0.00 11.36915    crim      0 116
    # 12   12      9        14         15  126.15810 0.00     7.57   lstat      0 212
    # 13   13      9        NA         NA   20.83679   NA   <leaf>     nox      1  34
    # 14   14     12        NA         NA   12.63760   NA   <leaf>     dis      1  42
    # 15   15     12        NA         NA   66.02809   NA   <leaf>    crim      1 170
    # 16   16     11        NA         NA   32.28858   NA   <leaf>   lstat      1  77
    # 17   17     11        NA         NA   76.00906 0.02   <leaf>     nox      1  39
    #        yval
    # 1  22.53281
    # 2  21.11071
    # 3  45.09667
    # 4  24.45924
    # 5  14.62037
    # 6  22.73242
    # 7  32.08103
    # 8  33.13000
    # 9  22.30976
    # 10 18.32826
    # 11 13.15000
    # 12 21.68821
    # 13 26.18529
    # 14 23.95000
    # 15 21.12941
    # 16 14.35195
    # 17 10.77692

## Classification Tree

    htt_iris <- HTT(X = iris[, 1:4], Y = iris[, 5], controls = htt_control(pt = 0.01))
    # prediction 
    table(predict(htt_iris), iris[, 5])

    #             
    #              setosa versicolor virginica
    #   setosa         50          0         0
    #   versicolor      0         49         5
    #   virginica       0          1        45

## Multivariate regression Tree

    data("ENB")
    set.seed(1)
    train = sample(1:nrow(ENB), floor(nrow(ENB)*0.8))
    train_x = ENB[train, 1:8]
    train_y = ENB[train, 9:10]
    test_x = ENB[-train, 1:8]
    test_y = ENB[-train, 9:10]
    htt_enb = HTT(train_x, train_y, controls = htt_control(pt = 0.05))
    # prediction
    pred = predict(htt_enb, newdata = test_x)
    # MAE
    colMeans(abs(pred - test_y))

    #        Y1        Y2 
    # 0.6244703 1.2953123

    # MSE
    colMeans(abs(pred - test_y)^2)

    #       Y1       Y2 
    # 1.239398 3.668416
