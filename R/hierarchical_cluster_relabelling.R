
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
#' Convert/test input to a standard format for cluster relabelling with
#'     \code{\link{relabel_hcluslabels}}.
#'
#' @param x input should be convertible to character matrix such that
#'     rows are labelled observations while columns contain hierarchical labels
#'     of increasing complexity.
#' @export
as_hcluslabels <- function(x) {
  UseMethod("as_hcluslabels")
}

#' @describeIn as_hcluslabels Check if object is in a valid hcluslabels format
is_hcluslabels <- function(x) {
  tryCatch({
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

validate_hcluslabels <- function(x) {

  cls <- class(x)

  values <- unclass(x)
  labels <- apply(values, 2, function(kk) sort(unique(kk)))

  label_lengths <- vapply(labels, FUN = length, FUN.VALUE = 1L)

  if ( nrow(values) < ncol(values) ) {
    stop("Error: cannot have more clusters (columns) than observations (rows)")
  }

  if ( ! "hcluslabels" %in% cls ) {
    stop("Error: not a hcluslabels object.")
  }

  if ( !all(diff(label_lengths) >= 1) ) {
    stop(
      "Error: number of unique levels must increase monotonically with column",
      call. = FALSE
    )
  }
  x
}


apply_lsap_labels <- function(x, lut) {
  unname(lut[match(x, names(lut))])
}

#' lsap_relabel_pair
#'
#' Function creates and applies a mapping from new -> old labels using
#'     a linear sum assignment solver (LSAP, see:
#'     \code{\link[clue]{solve_LSAP}}).
#'     For application to a set of hierarchical cluster labels see
#'     \code{\link{relabel_hcluslabels}}.
#'
#' @param new A vector of labels for subjects, to be modified such that labels
#'     in new are matched to best fitting labels in old.
#' @param old A vector of labels for subjects, values of old will be used to
#'     relabel new.
#' @return A copy of \code{new} with the updated labels.
#' @examples
#' lsap_relabel_pair(c("A", "B", "C"), c("x", "y", "z"))
#' @importFrom clue solve_LSAP
#' @export
lsap_relabel_pair <- function(new, old) {

  new <- as.character(new)
  old <- as.character(old)

  tab <- table(old, new)

  cn <- colnames(tab)
  rn <- rownames(tab)

  # solve_LSAP result is vector 1:nrow(tab) containing column matches:
  res <- clue::solve_LSAP(tab, maximum = TRUE)
  # extract matching colnames:
  res <- setNames(rn[match(seq_along(cn), res)], cn)

  # if no NA then we are working with k1==k2 and early return:
  if ( !any(is.na(res)) ) {
    return(apply_lsap_labels(new, res))
  }
  # otherwise we need to replace NA values with elements of cn not matched to
  #   elements of rn:
  possibles <- cn[cn %!in% rn]
  # Check we have enough:
  if ( length(possibles) < sum(is.na(res)) ) {
    stop("not enough names.")
  }
  # Assign the non-matched labels arbitrarily:
  res[is.na(res)] <- possibles[seq_along(res[is.na(res)])]
  # return:
  apply_lsap_labels(new, res)
}

#'  relabel_hcluslabels
#'
#'  Function to harmonise labels for clustering results on the same sample
#'      with increasing number of partitions - termed hierarchical clustering.
#'      For pairwise relabelling and non-hierarchical label matching
#'      (i.e. k1 == k2), see \code{\link{lsap_relabel_pair}}.
#'
#'  Clustering algorithms that accept a parameter k for the exact number of
#'  clusters to estimate will naturally return a variable number of clusters
#'  depending on that parameter value. Each increase in k produces a new label
#'  and some subjects will be assigned to the new label.
#'
#'  In addition to this, cluster labels are (commonly) assigned arbitrarily
#'  such that even when the exact same cluster structure is identified,
#'  which cluster is labelled as "A" vs. "B" is arbitrary. To illustrate, the
#'  following two clusterings of 10 patients into 4 clusters A,...,D is
#'  equivalent:
#'
#'
#'     A, A, A, B, C, C, C, C, D, D
#'
#'
#'     C, C, C, A, B, B, B, B, D, D
#'
#'  If we map the labels A -> B, B -> C, C -> A, and D -> D the two labellings
#'  are identical.
#'
#'  Given this setting, when evaluating results for a clustering algorithm at
#'  different values of k we don't want the assignment of labels to vary
#'  drastically from (e.g.) k=2 to k=3. Rather we want the cluster in k=3 which
#'  most resembles cluster "A" in k=2 to be labelled "A" (and not "B" or "C").
#'
#'  Note: this assumes a degree of continuity between
#'  clusterings at successive k values. It will not be useful if clusterings
#'  are at random. However in real data the two clusters found at k=2 rarely
#'  split at random into the 3 clusters at k=3 and it is more typical to have
#'  one of the k=2 labels split more into the new cluster.
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
#'  Without the relabelling the analyst needs to recall that "B" in k=3 was
#'  labelled "A" in the k=2 solution, which adds cognitive load to
#'  interpretation.
#'
#'  This function uses the \code{\link[clue]{solve_LSAP}} function from the
#'      clue package. See Hornik 2005 J. Stat. Software
#'      (\url{https://www.jstatsoft.org/article/view/v014i12})
#'
#' @param x matrix of clustering labels organised with subjects as rows,
#'        k-solutions as columns. Must be ordered with increasing complexity
#'        of solutions from left to right.
#' @return a matrix of class hcluslabels containing relabelled clusters from x
#' @importFrom stats setNames
#' @export
#' @examples
#' # Make some example data:
#' mat <- cbind(rep(c("A"), times = c(30)),
#'             rep(c("A","B"), times = c(10,20)),
#'             rep(c("C","B","A"), times = c(10,10,10)))
#' # Show that labels do not line up well:
#' #   most labelled A at k=1 are labelled B at k=2:
#' table(mat[,2], mat[,1], dnn = c("k=2", "k=1"))
#' #   the k=3 labels are jumbled:
#' table(mat[,3], mat[,2], dnn = c("k=3", "k=2"))
#'
#' # Apply the reordering:
#' rmat <- relabel_hcluslabels(mat)
#'
#' # inspect the remapping of labels:
#' attr(rmat, "mapping")
#'
#' # view the data:
#' rmat
#'
relabel_hcluslabels <- function(x) {
  x <- as_hcluslabels(x)

  x_orig <- x

  nk <- ncol(x)

  for ( i in seq(from = 2, to = nk) ) {
    x[, i] <- lsap_relabel_pair(x[, i], x[, i - 1])
  }

  mapping <- lapply(1:nk,
                    function(k) {
                      table(x[, k], x_orig[, k], dnn = c("new", "old"))
                    })
  names(mapping) <- colnames(x)

  return(structure(x, mapping = mapping))
}

#' check_relabelling
#'
#' Given a hcluslabels matrix check whether the relabelling was successful.
#'
#' @param new a matrix of recoded labels, the result of relabel_hcluslabels(old)
#' @param old an optional matrix of raw labels
#' @return unchanged labels
#' @export
check_relabelling <- function(new, old = NULL) {

  if ( ! "mapping" %in% names(attributes(new)) ) {
    stop("new does not have the expected 'mapping' attribute.")
  }

  map <- attr(new, "mapping")
  # Check that each mapping is unique.
  #   a mapping is unique if there is one non-zero entry for every row and
  #   every column.
  map_chk <- lapply(map, function(x) x > 0)

  chk_unique_mapping <- function(a) {
    nc <- ncol(a)
    nr <- nrow(a)
    chk_r <- all.equal(as.vector(rowSums(a)),
                       rep(1, nc))
    chk_c <- all.equal(as.vector(colSums(a)),
                       rep(1, nr))
    !chk_r | !chk_c
  }

  lapply(map_chk,
         function(x) {
           if ( chk_unique_mapping(x) ) {
             stop("mapping margins wrong")
           }
         } )

  if ( !is.null(old) ) {
    if ( !identical(dim(new), dim(old)) ) {
      stop("Dimensions do not match")
    }
  }
  return(new)
}
