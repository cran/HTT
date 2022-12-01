print.htt <- function(x, ...) {
  if (!inherits(x, "htt"))
    stop("Not a legitimate \"htt\" object")
  frame <- x$frame
  var.type <- x$var.type
  var.name <- names(var.type)
  cat("     Hypothesis Testing Tree \n\n")
  cat("node, split, n, pvalue\n")
  cat("* denotes terminal node\n\n")
  seq <- travel(frame, 1)
  cat("[", 1, "] ", "root", sep = "")
  cat(paste0("   (n = ", frame[1, "n"], ",",
             " pvalue = ", frame[1, "pval"], ")\n", seq = ""))
  for (num in seq[-1]) {
    cat(paste0(rep("|  ", times = BTHeight(frame, num) - 1)),
        "[", num, "] ", sep = "")
    parent <- frame[num, "parent"]
    if (!var.type[frame[parent, "var"]]) {
      if (num == frame[parent, "leftChild"]) {
        cat(frame[parent, "var"], "<=", frame[parent, "split"], sep = "")
      } else {
        cat(frame[parent, "var"], ">", frame[parent, "split"], sep = "")
      }
    } else {
      sp <- paste0(strsplit(frame[parent, "split"], "}")[[1]], "}")
      if (num == frame[parent, "leftChild"]) {
        cat(frame[parent, "var"], " == ", sp[1], sep = "")
      } else {
        cat(frame[parent, "var"], " == ", sp[2], sep = "")
      }
    }
    if (x$frame[num, "isleaf"] == 0) {
      cat(paste0("   (n = ", frame[num, "n"], ",", " pvalue = ",
                 x$frame[num, "pval"], ")\n", seq = ""))
    } else {
      cat(paste0("   (n = ", frame[num, "n"], ")", " *", "\n", seq = ""))
    }
  }
  invisible(x)
}
