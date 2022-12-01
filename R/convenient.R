### find all leafs below the node
node_all <- function(frame, node) {
  if (frame[node, "isleaf"] == 1) {
    return(node)
  } else {
    left <- node_all(frame = frame, node = frame[node, "leftChild"])
    right <- node_all(frame = frame, node = frame[node, "rightChild"])
    return(c(left, right))
  }
}
### height of a node in the tree
tree.height <- function(frame, node) {
  if (frame[node, "parent"] == 0) {
    return(1)
  } else {
    parentHeight <- tree.height(frame, node = frame[node, "parent"])
    return(parentHeight + 1)
  }
}
travel <- function(frame, node) {
  if (frame[node, "isleaf"] == 1) {
    return(node)
  } else {
    left <- travel(frame = frame, node = frame[node, "leftChild"])
    right <- travel(frame = frame, node = frame[node, "rightChild"])
    return(c(node, left, right))
  }
}
BTHeight <- function(frame, node) {
  if (frame[node, "parent"] == 0) {
    return(1)
  } else {
    parentHeight <- BTHeight(frame, node = frame[node, "parent"])
    return(parentHeight + 1)
  }
}
