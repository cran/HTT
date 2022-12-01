printsplit <- function(object) {
  if (!inherits(object, "htt"))
    stop("Not a legitimate \"htt\" object")
  ff <- object$frame
  leaf <- which(ff$isleaf == 1)
  ff[leaf, c("leftChild", "rightChild", "pval")] <- "*"
  ff[is.na(ff$var), "var"] <- "*"
  print(ff)
}
