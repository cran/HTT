predict.htt <- function(object, newdata,
                        type = c("response", "prob", "node"), ...)
{
  if (!inherits(object, "htt"))
    stop("Not a legitimate \"htt\" object")
  if (missing(newdata)) {
    newdata <- object$X
  }
  if (!is.data.frame(newdata)) {
    stop("newdata should be a dataframe, not matrix or array.")
  }
  frame <- object$frame
  used <- sort(unique(na.omit(frame$var)))
  name.new <- colnames(newdata)
  ind <- which(used %in% name.new == 0)
  if (length(ind) > 0) {
    cat(paste(used[ind], collapse = ", "), "are not in newdata")
    stop("Incomplete newdata.")
  }
  if (missing(type)) {
    type <- "response"
  }
  if (type == "prob") {
    if (object$method != "classification") {
      Stop("prob type is not sitable for regression tree.")
    }
  }
  X <- object$X
  var.type <- object$var.type
  Stop <- 0
  for (var in used) {
    if (var.type[var] == 0) {
      if (!is.numeric(newdata[, var])) {
        cat(var, "should be continuous.\n")
        Stop <- 1
      }
    } else if (var.type[var] == 1) {
      if (!is.ordered(newdata[, var])) {
        cat(var, "should be ordinal\n")
        Stop <- 1
      }
      set_diff <- setdiff(levels(newdata[, var]), levels(X[, var]))
      if (length(set_diff) > 0) {
        cat("factor", var, "has new levels",
            paste(set_diff, collapse = ", "), "\n")
        Stop <- 1
      }
    } else {
      if (!is.factor(newdata[, var]) | is.ordered(newdata[, var])) {
        cat(var, "should be nominal\n")
        Stop <- 1
      }
      set_diff <- setdiff(levels(newdata[, var]), levels(X[, var]))
      if (length(set_diff) > 0) {
        cat("factor", var, "has new levels",
            paste(set_diff, collapse = ", "), "\n")
        Stop <- 1
      }
    }
  }
  if (Stop == 1) {
    stop("Incorrect newdata.")
  }
  pred_label <- rep(1, nrow(newdata))
  index <- 1:nrow(newdata)
  while (length(index) > 0) {
    uni <- unique(pred_label[index])
    for (i in uni) {
      ind <- index[which(pred_label[index] == i)]
      if (frame[i, "isleaf"] == 1) {
        index <- setdiff(index, ind)
      } else {
        split <- frame[i, "split"]
        if (var.type[frame[i, "var"]] == 0) {
          split <- as.numeric(split)
          ind_left <- ind[which(newdata[ind, frame[i, "var"]] <= split)]
          pred_label[ind_left] <- frame[i, "leftChild"]
          ind_right <- ind[which(newdata[ind, frame[i, "var"]] > split)]
          pred_label[ind_right] <- frame[i, "rightChild"]
        } else {
          split <- strsplit(split, "}")[[1]]
          rsplit <- strsplit(gsub("\\{", "", split[2]), ",")[[1]]
          ind_right <- ind[which(newdata[ind, frame[i, "var"]]%in%rsplit)]
          pred_label[ind_right] <- frame[i, "rightChild"]
          ind_left <- setdiff(ind, ind_right)
          pred_label[ind_left] <- frame[i, "leftChild"]
        }
      }
    }
  }
  if (type == "response") {
    yval <- frame[, grep("yval", names(frame))]
    yval <- as.matrix(yval)
    pred_val <- yval[pred_label, ]
    return(response = pred_val)
  } else if (type == "node") {
    return(node = pred_label)
  } else {
    yprob <- frame[, grep("prob", names(frame))]
    prob <- yprob[pred_label, ]
    return(prob = prob)
  }
}
