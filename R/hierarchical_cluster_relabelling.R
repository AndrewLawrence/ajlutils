# Rationale:


new_hcluslabels <- function(x) {
  stopifnot(is.matrix(x))
  if ( !is.character(x) ) {
    x <- apply(x, 2, as.character)
  }
  structure(x,
            class = c("hcluslabels", "matrix", "array"))
}

#' as_hcluslabels
#'
#' Convert/test input standard format for cluster relabeling.
#'
#' @param x input should be convertable to character matrix
#' @export
as_hcluslabels <- function(x) {
  UseMethod("as_hcluslabels")
}

#' @describeIn as_hcluslabels Check if object is in a valid hcluslabels format
is_hcluslabels <- function(x) {
  tryCatch( {
    as_hcluslabels(x)
    TRUE
    },
    error = function(e) FALSE)
}


#' @export
as_hcluslabels.default <- function(x) {
  validate_hcluslabels(new_hcluslabels(as.matrix(x)))
}

#' @export
as_hcluslabels.matrix <- function(x) {
  validate_hcluslabels(new_hcluslabels(x))
}

#' @export
as_hcluslabels.data.frame <- function(x) {
  as_hcluslabels.matrix(as.matrix(x))
}

#' @export
as_hcluslabels.hcluslabels <- function(x) {
  validate_hcluslabels(x)
}

get_labels <- function(x) {
  apply(x, 2, function(kk) sort(unique(kk)))
}

validate_hcluslabels <- function(x) {

  cls <- class(x)

  values <- unclass(x)
  labels <- get_labels(values)

  label_lengths <- vapply(labels, FUN = length, FUN.VALUE = 1L)

  if ( nrow(values) < ncol(values) ) {
    stop("Error: cannot have more clusters (columns) than observations (rows)")
  }

  if ( ! "hcluslabels" %in% cls ) {
    stop("Error: not a hcluslabels object. Try ")
  }

  if ( !all(diff(label_lengths) >= 1) ) {
    stop(
      "Error: number of unique levels must increase monotonically with column",
      call. = FALSE
    )
  }
  x
}

cluster_relabel <- function(labs2, labs1) {
  # note: strictly labs2 must have more labels than labs1

  # joint tabulation:
  xt <- table(labs1, labs2)
  # express as a proportion:
  p_xt <- prop.table(xt, margin = 1)

  # for each row (i.e. previous label), find the column with the
  #   greatest proportion overlap:
  max_idx <- apply(xt, 1, which.max) # index of max
  # turn this into a lookuptable:
  lut <- setNames(colnames(p_xt)[unname(max_idx)],
                      rownames(p_xt))

  # apply the lut:
  newlabs <- names(lut)[match(colnames(p_xt), unname(lut))]

  # fill in NA (keep new labels):
  newlabs <- setNames(
    ifelse(is.na(newlabs), colnames(p_xt), newlabs),
    colnames(p_xt))

  # apply this to the data:
  return(unname(newlabs[match(labs2, names(newlabs))]))
  return(list(p_xt = p_xt,
              max_idx = max_idx,
              lut = lut,
              newlabs = newlabs))
}


#'  relabel_hcluslabels
#'
#'  Function to relabel hierarchical clustering results
#'
#'  Clustering algorithms that accept a parameter k for the exact number of
#'  clusters to estimate will naturally return a variable number of clusters
#'  depending on the parameter value.
#'
#'  Clusters labels are typically assigned arbitrarily by clustering algorithms
#'  To illustrate, assigning 4 clusters to 10 patients in the following ways
#'  is equivalent as we can map the labels A -> B, B -> C, C -> A, and D -> D to
#'  obtain the same result:
#'
#'
#'     A, A, A, B, C, C, C, C, D, D
#'
#'
#'     C, C, C, A, B, B, B, B, D, D
#'
#'  Given this setting, when evaluating results for a clustering algorithm at
#'  different values of k we don't want the labels to vary drastically
#'  from (e.g.) k=2 to k=3 rather we want the cluster in k=3 which
#'  most resembles cluster "A" in k=2 to be labelled "A" and not "B" or "C".
#'  But this behaviour is not possible given only the input to the clustering
#'  algorithm. Note: this is only useful if there is continuity between
#'  clusterings at successive k values. It will not be useful if clusterings
#'  are at random. In real data the two clusters found at k=2 rarely split
#'  at random into the 3 clusters at k=3.
#'
#'  To illustrate, in the scenario:
#'
#'   k=2: A A A B B B B B B
#'
#'   k=3: B B B A A A A C C
#'
#'   we would prefer to adjust k=3 labels such that:
#'
#'   k=2: A A A B B B B B B
#'
#'   k=3: A A A B B B B C C
#'
#'  The relabelling makes it easier to see that:
#'
#'     1) The "A" class is preserved from k=2 solution to the k=3 solution
#'
#'     2) The "B" class splits into two clusters "B" and "C"
#'
#'  Without the relabelling the analyst must recall that "B" in k=3 was labelled
#'  "A" in the k=2 solution, which adds cognitive load to interpretation.
#'
#'  This function implements a simple greedy algorithm to relabel clusters in a
#'  cascading manner from left to right.
#'
#' @param x matrix of clustering labels organised with subjects as rows,
#'        k-solutions as columns with increasing complexity solutions
#'        from left to right.
#' @importFrom stats setNames
#' @export
relabel_hcluslabels <- function(x) {
  x <- as_hcluslabels(x)

  nk <- ncol(x)

  for ( i in seq(from = 2, to = nk) ) {
    x[, i] <- cluster_relabel(x[, i], x[, i - 1])
  }
  x
}

# HERE: bug in relabel_hcluslabels function:
tmp <- structure(c("Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z",
                   "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "Z", "S",
                   "L", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "L",
                   "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "L",
                   "L", "L", "S", "S", "S", "S", "S", "L", "S", "S", "S", "S", "S",
                   "S", "L", "L", "S", "S", "S", "S", "S", "L", "L", "L", "S", "S",
                   "S", "L", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S",
                   "S", "L", "S", "S", "S", "L", "S", "S", "S", "S", "S", "S", "S",
                   "L", "S", "S", "L", "S", "L", "S", "L", "S", "S", "S", "S", "S",
                   "L", "S", "S", "L", "S", "S", "S", "S", "X", "H", "X", "X", "X",
                   "X", "L", "X", "X", "X", "X", "L", "L", "H", "X", "X", "L", "X",
                   "X", "X", "L", "L", "X", "X", "X", "X", "H", "H", "H", "L", "X",
                   "L", "L", "X", "H", "X", "X", "L", "X", "X", "X", "H", "H", "X",
                   "X", "X", "L", "X", "H", "H", "H", "L", "X", "X", "H", "X", "L",
                   "X", "L", "X", "L", "X", "X", "X", "X", "L", "X", "H", "X", "X",
                   "X", "H", "X", "L", "L", "X", "X", "X", "L", "H", "X", "X", "H",
                   "X", "H", "X", "H", "X", "L", "X", "X", "X", "H", "X", "X", "H",
                   "X", "X", "X", "X", "D", "W", "D", "D", "D", "D", "H", "A", "D",
                   "D", "D", "H", "H", "W", "D", "A", "H", "D", "D", "D", "H", "H",
                   "A", "D", "D", "D", "W", "W", "W", "H", "A", "H", "H", "D", "W",
                   "A", "D", "H", "D", "A", "D", "W", "W", "D", "D", "A", "H", "D",
                   "W", "W", "W", "H", "A", "A", "W", "D", "H", "D", "H", "D", "H",
                   "D", "D", "A", "A", "H", "D", "W", "D", "A", "D", "W", "D", "H",
                   "H", "D", "D", "D", "H", "W", "A", "A", "W", "D", "W", "D", "W",
                   "A", "H", "A", "A", "A", "W", "D", "D", "W", "D", "A", "A", "D",
                   "L", "N", "L", "S", "L", "S", "E", "Z", "L", "L", "S", "E", "E",
                   "N", "S", "Z", "E", "S", "S", "L", "E", "E", "Z", "L", "L", "S",
                   "N", "N", "N", "E", "Z", "E", "E", "L", "N", "Z", "S", "E", "S",
                   "Z", "S", "N", "N", "L", "L", "Z", "E", "S", "N", "N", "N", "E",
                   "Z", "Z", "N", "L", "E", "S", "E", "S", "E", "L", "L", "Z", "Z",
                   "E", "L", "N", "S", "Z", "L", "N", "S", "E", "E", "L", "L", "L",
                   "E", "N", "Z", "Z", "N", "S", "N", "L", "N", "Z", "E", "Z", "Z",
                   "Z", "N", "S", "S", "N", "S", "Z", "Z", "S"), dim = c(100L, 5L
                   ))

apply(relabel_hcluslabels(tmp), 2, table)
# [[1]]
#
#   Z
# 100
#
# [[2]]
#
#  L  Z
# 20 80
#
# [[3]]
#
#  L  Z
# 40 60
#
# [[4]]
#
#  A  L  W  Z
# 20 20 20 40
#
# [[5]]
#
#  A  L  S  W  Z
# 20 20 20 20 20
#
