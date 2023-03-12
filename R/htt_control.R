### control the hyper parameters
htt_control <- function(teststat = c("energy0", "energy1"),
                        testtype = c("permutation", "fastpermutation"),
                        alpha = 1, pt = 0.05, minsplit = 30,
                        minbucket = round(minsplit/3),
                        R = 199, nmin = 1000) {
  teststat <- match.arg(teststat)
  testtype <- match.arg(testtype)
  RET <- list()
  if (teststat %in% teststat) {
    RET$teststat <- teststat
  } else {
    stop(sQuote("teststat"), teststat, " not defined")
  }
  if (testtype %in% testtype) {
    RET$testtype <- testtype
  } else {
    stop(sQuote("testtype"), testtype, " not defined")
  }
  if (0 < alpha && alpha <= 2) {
    RET$alpha <- alpha
  } else {
    stop(sQuote("alpha"), " shoule between 0 and 2")
  }
  if (0 <= pt & pt <= 1) {
    RET$pt <- pt
  } else {
    stop(sQuote("pt"), " shoule between 0 and 1")
  }
  if (minsplit >= 20) {
    RET$minsplit <- as.integer(minsplit)
  } else {
    stop(sQuote("minsplit"), " should be at least 20")
  }
  if (minbucket >= 7) {
    RET$minbucket <- as.integer(minbucket)
  } else {
    stop(sQuote("minbucket"), " should be at least 7")
  }
  RET$minsplit = max(RET$minsplit, 3*RET$minbucket)
  if (R >= 19) {
    RET$R <- as.integer(R)
  } else {
    stop(sQuote("R"), " should larger than 99")
  }
  if (nmin >= 100){
    RET$nmin <- as.integer(nmin)
  } else {
    stop(sQuote("nmin"), " should larger than 100")
  }
  RET
}
