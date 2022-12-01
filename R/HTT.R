HTT <- function(X, Y, method, distance = NULL, controls = htt_control()) {
  if (missing(method)) {
    if (is.factor(Y) || is.character(Y)) {
      method <- "classification"
    } else if ((is.matrix(Y) || is.vector(Y)) & is.numeric(Y)) {
      method <- "regression"
    } else if (is.data.frame(Y) & (sum(sapply(Y, is.numeric)) == ncol(Y))) {
      Y = as.matrix(Y)
      method <- "regression"
    } else {
      stop("Y should be a numeric or factor!")
    }
  }
  method.int <- pmatch(method, c("regression", "classification"))
  if (is.na(method.int))
    stop("Invalid method")
  if (is.matrix(X) & is.numeric(X)) {
    XX <- X
    X <- data.frame(X)
    var_type <- rep(0, ncol(X))
    names(var_type) <- colnames(X)
  } else if (is.data.frame(X)) {
    var_type <- sapply(X, function(x) {
      if (is.numeric(x)) {
        return(0)
      } else if (is.ordered(x)) {
        return(1)
      } else if (is.factor(x)) {
        return(2)
      } else {
        stop(" 'X' should be a numeric vector/matrix, factor or dataframe.")
      }
    })
    XX <- lapply(X, function(x) {
      if (!is.factor(x)) {
        x
      } else {
        as.numeric(x)
      }
    })
    XX <- as.matrix(as.data.frame(XX))
  } else if (is.numeric(X)) {
    X <- data.frame(X)
    var_type <- 0
    names(var_type) <- colnames(X)
    XX <- as.matrix(X)
  } else if (is.factor(X)) {
    var_type <- 2
    XX <- as.matrix(as.numeric(X))
    X <- data.frame(X)
    names(var_type) <- colnames(X)
  } else {
    stop(" 'X' should be a numeric vector/matrix, factor or dataframe.")
  }
  if (nrow(X) != nrow(as.matrix(Y))) {
    stop("variable lengths differ (found for 'X' and 'Y')")
  }
  alpha <- controls$alpha
  controls$maxnode <- ceiling(nrow(XX)/(controls$minsplit)) * 4 + 1

  if (is.null(distance)) {
    if (method.int == 1L) {
      dmat <- distance(as.matrix(Y), alpha)
    }
    if (method.int == 2L) {
      dmat <- outer(Y, Y, FUN = function(x, y) x != y)
    }
    fit <- TreeGrow(XX, dmat, var_type, controls)
  } else {
    if (!is.matrix(distance)) {
      stop(" 'distance' should be a matrix")
    }
    if (nrow(distance) != ncol(distance)) {
      stop(" 'distance' should be a square matrix")
    }
    if (nrow(distance) != nrow(X)) {
      stop(" 'distance' should have the same rows as X")
    }
    if (!is.numeric(distance)) {
      stop(" 'distance' should be a numeric matrix")
    }
    if (any(distance < 0)) {
      stop(" 'distance' should be a non-negative matrix")
    }
    if (sum(diag(distance)) > 0) {
      stop(" The diagonal of 'distance' should be 0")
    }
    if (any(distance != t(distance))) {
      stop(" 'distance' should be a symmetric matrix")
    }
    fit <- TreeGrow(XX, distance, var_type, controls)
  }
  frame <- fit$frame
  where <- fit$where + 1
  if (!is.null(rownames(X))) {
    names(where) <- rownames(X)
  } else {
    names(where) <- 1:nrow(X)
  }
  n0 <- which.max(frame$node)
  frame <- frame[1:n0, ]
  frame[, c("node", "parent", "leftChild", "rightChild")] <-
    frame[, c("node", "parent", "leftChild", "rightChild")] + 1
  frame[1, "parent"] <- 0
  leaf <- which(frame$isleaf == 1)
  frame$leftChild[leaf] <- NA
  frame$rightChild[leaf] <- NA
  frame$statistic <- round(frame$statistic, getOption("digits") - 2)
  frame$split <- as.character(round(frame$split, getOption("digits")))
  frame$split[leaf] <- "<leaf>"
  ind <- which(is.na(frame$split) & frame$isleaf == 0)
  if (length(ind) > 0) {
    for (i in ind) {
      lchild <- frame[i, "leftChild"]
      rchild <- frame[i, "rightChild"]
      d <- frame[i, "var"]
      ind_l <- which(where %in% node_all(frame, lchild))
      ind_r <- which(where %in% node_all(frame, rchild))
      split_l <- sort(unique(X[ind_l, d]))
      split_r <- sort(unique(X[ind_r, d]))
      string <- paste0("{", paste0(split_l, collapse = ","), "}{",
                       paste0(split_r, collapse = ","), "}")
      frame[i, "split"] <- string
    }
  }
  if (!is.factor(Y) & is.vector(Y)) {
    yval <- sapply(1:nrow(frame), function(i) {
      Node <- node_all(frame, i)
      ind <- which(where %in% Node)
      return(mean(Y[ind]))
    })
    frame$yval <- round(yval, getOption("digits") - 2)
  } else if (!is.factor(Y) & is.matrix(Y)) {
    yval <- sapply(1:nrow(frame), function(i) {
      Node <- node_all(frame, i)
      ind <- which(where %in% Node)
      return(colMeans(Y[ind, ]))
    })
    yval <- t(yval)
    colnames(yval) <- paste0("yval", 1:ncol(yval))
    frame <- data.frame(frame, round(yval, getOption("digits") - 2))
  } else {
    y <- droplevels(Y)
    k <- nlevels(y)
    Prob <- matrix(0, nrow(frame), k)
    colnames(Prob) <- paste0("prob_", levels(y))
    yval <- rep("", nrow(frame))
    for (i in 1:nrow(frame)) {
      Node <- node_all(frame, i)
      ind <- which(where %in% Node)
      n <- length(ind)
      prob <- table(y[ind])/n
      Prob[i, ] <- prob
      yval[i] <- names(which.max(prob))
    }
    frame <- data.frame(frame, yval, round(Prob, 3))
  }
  frame$var <- colnames(X)[frame$var]
  res <- list(frame = frame, where = where, method = method, control = controls,
              X = X, var.type = var_type)
  class(res) <- "htt"
  invisible(res)
  return(res)
}
