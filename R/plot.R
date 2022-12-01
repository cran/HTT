plot.htt <- function(x, digits = 3, line.color = "blue", node.color = "black",
                     line.type = c("straight", "curved"),
                     layout = c("tree", "dendrogram"), ...) {
  if (!inherits(x, "htt"))
    stop("Not a legitimate \"htt\" object")
  line.type <- match.arg(line.type)
  layout <- match.arg(layout)
  frame <- x$frame
  var.type <- x$var.type
  n <- nrow(frame)
  yval <- frame[, grep("yval", names(frame))]
  yval <- as.matrix(yval)
  if (ncol(yval) > 1) {
    yval <- round(yval, digits)
    p <- ncol(yval)
    yval <- apply(yval, 1, function(x)
      paste0("y", 1:p, "=", x, collapse = "\n"))
  } else {
    yval <- paste0("y=", c(yval))
  }
  node_info <- paste0("node", 1:n, "\n", yval)
  terminal <- which(frame$isleaf == 1)
  pval <- round(frame$pval[-terminal], digits)
  pval <- ifelse(pval == 0, paste0("p < ", 1/(1 + x$control$R)),
                 paste0("p = ", pval))
  node_info[-terminal] <- paste0(node_info[-terminal], "\n", pval)
  node_info[terminal] <- paste0(node_info[terminal], "\n", "n = ",
                                frame$n[terminal])
  if (n == 1) {
    from <- 1
    to <- 1
    label.name <- ""
  } else {
    interior <- which(frame$isleaf == 0)
    from <- rep(interior, each = 2)
    to_L <- frame$leftChild[interior]
    to_R <- frame$rightChild[interior]
    to <- rep(0, length(from))
    to[seq(1, length(to), by = 2)] <- to_L
    to[seq(2, length(to), by = 2)] <- to_R
    # to = c(to_L, to_R)
    sp_L <- vector("character", length(interior))
    sp_R <- vector("character", length(interior))
    for (i in 1:length(interior)) {
      var <- frame[interior[i], "var"]
      split <- frame[interior[i], "split"]
      if (!var.type[var]) {
        split <- round(as.numeric(split), digits)
        sp_L[i] <- paste0(var, "<=", split)
        sp_R[i] <- paste0(var, ">", split)
      } else {
        sp <- paste0(strsplit(split, "}")[[1]], "}")
        sp_L[i] <- paste0(var, "==", sp[1])
        sp_R[i] <- paste0(var, "==", sp[2])
      }
    }
    sp_L <- paste0(sp_L)
    sp_R <- paste0(sp_R, "\n")
    label.name <- rep("", length(to))
    label.name[seq(1, length(label.name), by = 2)] <- sp_L
    label.name[seq(2, length(label.name), by = 2)] <- sp_R
    # label.name = c(sp_L, sp_R)
  }
  edge <- cbind(from = node_info[from],
                to = node_info[to],
                label.name = label.name)
  gr <- graph_from_data_frame(edge)
  # name = node_info
  name <- 1
  if (line.type == "straight") {
    # set_graph_style(plot_margin = margin(1,1,1,1))
    ggraph(gr, layout = layout, circular = F) +
      geom_node_point() +
      geom_edge_link(aes(label = label.name), n = 10, color = line.color) +
      geom_node_label(aes(label = name), color = node.color) + theme_void()
  } else {
    # set_graph_style(plot_margin = margin(1,1,1,1))
    ggraph(gr, layout = layout, circular = F) +
      geom_node_point() +
      geom_edge_diagonal(aes(label = label.name), color = line.color) +
      geom_node_label(aes(label = name), color = node.color) + theme_void()
  }
}
