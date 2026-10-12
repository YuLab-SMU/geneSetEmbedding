#' Clean a PPI edge list before building a graph
#'
#' Optional preprocessing for edge tables: drop self-loops, collapse undirected
#' duplicates to one orientation, and optionally keep edges above a weight
#' cutoff. Call this before \code{\link{gsemb_build_graph}} /
#' \code{\link{gsemb_fit}} when the input may contain self-loops, both
#' directions of an undirected pair, or low-confidence scores. Already-clean
#' edge lists can skip this step. Rows with missing or empty endpoint IDs are
#' skipped with a warning giving the number of skipped rows.
#'
#' @param edges A data.frame containing at least two columns for endpoints.
#' @param node1,node2 Column names in \code{edges} for source/target.
#' @param weight Optional column name in \code{edges} for edge weights. Required
#'   when \code{score_cutoff} is not \code{NULL}.
#' @param score_cutoff Optional lower bound on \code{weight}. Rows with
#'   \code{weight <= score_cutoff} are dropped. Default \code{NULL} keeps all
#'   weights. For STRING \code{combined_score}, a common choice is \code{700}.
#' @param drop_self_loops Logical; drop rows where the two endpoints are equal.
#' @param undirected_unique Logical; for each undirected pair, keep a single
#'   row with endpoints ordered so \code{node1 < node2}. If both directions (or
#'   other duplicates) remain after that, keep the row with the largest
#'   \code{weight} when \code{weight} is set, otherwise the first row.
#'
#' @return A data.frame with the same columns as \code{edges}, possibly fewer
#'   rows. Endpoint columns are character.
#' @examples
#' edges <- data.frame(
#'   node1 = c("B", "A", "A", "C"),
#'   node2 = c("A", "B", "A", "D"),
#'   weight = c(800, 800, 900, 500)
#' )
#' gsemb_clean_edges(edges, weight = "weight", score_cutoff = 700)
#' @seealso \code{\link{gsemb_build_graph}}, \code{\link{gsemb_fit}}
#' @export
gsemb_clean_edges <- function(edges,
                              node1 = "node1",
                              node2 = "node2",
                              weight = NULL,
                              score_cutoff = NULL,
                              drop_self_loops = TRUE,
                              undirected_unique = TRUE) {
  if (!is.data.frame(edges)) stop("edges must be a data.frame")
  if (!node1 %in% names(edges)) stop("node1 column not found")
  if (!node2 %in% names(edges)) stop("node2 column not found")
  if (!is.null(score_cutoff) && is.null(weight)) {
    stop("weight must be set when score_cutoff is not NULL")
  }
  if (!is.null(weight) && !weight %in% names(edges)) stop("weight column not found")

  out <- edges
  a <- as.character(out[[node1]])
  b <- as.character(out[[node2]])
  out[[node1]] <- a
  out[[node2]] <- b

  missing_endpoint <- is.na(a) | is.na(b) | a == "" | b == ""
  n_missing <- sum(missing_endpoint)
  if (n_missing > 0L) {
    warning(sprintf("Dropped %d edges with missing endpoints", n_missing), call. = FALSE)
    out <- out[!missing_endpoint, , drop = FALSE]
    a <- out[[node1]]
    b <- out[[node2]]
  }

  if (isTRUE(drop_self_loops)) {
    out <- out[a != b, , drop = FALSE]
    a <- out[[node1]]
    b <- out[[node2]]
  }

  if (!is.null(score_cutoff)) {
    w <- as.numeric(out[[weight]])
    out <- out[is.finite(w) & w > score_cutoff, , drop = FALSE]
    a <- out[[node1]]
    b <- out[[node2]]
  }

  if (isTRUE(undirected_unique) && nrow(out) > 0) {
    left <- pmin(a, b)
    right <- pmax(a, b)
    out[[node1]] <- left
    out[[node2]] <- right
    pair <- paste(left, right, sep = "\r")
    if (!is.null(weight)) {
      w <- as.numeric(out[[weight]])
      w[!is.finite(w)] <- -Inf
      # Keep the highest-weight row per undirected pair.
      ord <- order(pair, -w, seq_len(nrow(out)))
      out <- out[ord, , drop = FALSE]
      pair <- pair[ord]
    }
    out <- out[!duplicated(pair), , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}

#' Build a PPI graph adjacency matrix
#'
#' Construct a sparse adjacency matrix from an edge list. Node names are stored
#' in the matrix dimnames and are used throughout the package. This function
#' does not drop self-loops, collapse bidirectional undirected pairs, or apply
#' a score cutoff; use \code{\link{gsemb_clean_edges}} first when needed.
#'
#' @param edges A data.frame containing at least two columns for endpoints.
#' @param node1,node2 Column names in \code{edges} for source/target.
#' @param weight Optional column name in \code{edges} providing edge weights.
#' @param directed Logical; whether to keep the graph directed.
#' @param nodes Optional character vector of node IDs to keep/order.
#'
#' @return A sparse adjacency matrix of class \code{dgCMatrix}.
#' @examples
#' # Create a small edge list
#' edges <- data.frame(
#'   node1 = c("A", "B", "C", "A"),
#'   node2 = c("B", "C", "A", "C"),
#'   weight = c(1.0, 2.0, 0.5, 1.5)
#' )
#'
#' # Build undirected weighted graph
#' adj <- gsemb_build_graph(edges, weight = "weight")
#' adj
#'
#' # Build directed unweighted graph
#' adj_dir <- gsemb_build_graph(edges, directed = TRUE)
#' adj_dir
#' @seealso \code{\link{gsemb_clean_edges}}
#' @export
gsemb_build_graph <- function(edges,
                              node1 = "node1",
                              node2 = "node2",
                              weight = NULL,
                              directed = FALSE,
                              nodes = NULL) {
  if (!is.data.frame(edges)) stop("edges must be a data.frame")
  if (!node1 %in% names(edges)) stop("node1 column not found")
  if (!node2 %in% names(edges)) stop("node2 column not found")

  src <- as.character(edges[[node1]])
  dst <- as.character(edges[[node2]])
  if (is.null(weight)) {
    w <- rep(1, length(src))
  } else {
    if (!weight %in% names(edges)) stop("weight column not found")
    w <- as.numeric(edges[[weight]])
    w[is.na(w)] <- 0
  }

  if (is.null(nodes)) {
    nodes <- sort(unique(c(src, dst)))
  } else {
    nodes <- unique(as.character(nodes))
  }
  idx1 <- match(src, nodes)
  idx2 <- match(dst, nodes)
  ok <- !is.na(idx1) & !is.na(idx2) & w != 0
  idx1 <- idx1[ok]
  idx2 <- idx2[ok]
  w <- w[ok]

  adj <- Matrix::sparseMatrix(i = idx1, j = idx2, x = w, dims = c(length(nodes), length(nodes)))
  if (!directed) {
    adj <- adj + Matrix::t(adj)
  }
  dimnames(adj) <- list(nodes, nodes)
  Matrix::drop0(adj)
}

#' Construct a random-walk transition matrix
#'
#' Normalize an adjacency matrix into a transition matrix by column- or
#' row-normalization.
#'
#' @param adj A sparse/dense matrix (typically produced by \code{gsemb_build_graph}).
#' @param normalize One of \code{"col"} (column-stochastic) or \code{"row"} (row-stochastic).
#' @param eps Small constant to avoid division by zero.
#'
#' @return A sparse transition matrix.
#' @examples
#' # Build a small adjacency matrix
#' edges <- data.frame(
#'   node1 = c("A", "B", "C"),
#'   node2 = c("B", "C", "A")
#' )
#' adj <- gsemb_build_graph(edges)
#'
#' # Column-normalized transition matrix
#' Wcol <- gsemb_transition_matrix(adj, normalize = "col")
#' Matrix::colSums(Wcol)
#'
#' # Row-normalized transition matrix
#' Wrow <- gsemb_transition_matrix(adj, normalize = "row")
#' Matrix::rowSums(Wrow)
#' @export
gsemb_transition_matrix <- function(adj, normalize = c("col", "row"), eps = 1e-12) {
  normalize <- match.arg(normalize)
  if (!inherits(adj, "Matrix")) stop("adj must be a Matrix sparse/dense matrix")
  if (normalize == "col") {
    s <- Matrix::colSums(adj)
    s <- pmax(as.numeric(s), eps)
    Dinv <- Matrix::Diagonal(x = 1 / s)
    W <- adj %*% Dinv
  } else {
    s <- Matrix::rowSums(adj)
    s <- pmax(as.numeric(s), eps)
    Dinv <- Matrix::Diagonal(x = 1 / s)
    W <- Dinv %*% adj
  }
  Matrix::drop0(W)
}

#' Select landmark nodes for diffusion features
#'
#' Choose landmark nodes by node strength, uniform random sampling,
#' weight-aware betweenness, or weight-aware K-Medoids. Edge cost is
#' \code{1 / weight}, so higher edge weight is a shorter path. For betweenness
#' and K-Medoids, nonzero edge weights must be finite and positive; the scale
#' may be 0--1, STRING 0--1000, or another positive scale.
#'
#' @param adj Adjacency matrix with rownames as node IDs.
#' @param k Number of landmarks to select (capped at number of nodes).
#' @param method Landmark selection method: \code{"degree"} ranks by the sum of
#'   incident edge weights, \code{"random"} samples nodes uniformly,
#'   \code{"betweenness"} ranks by betweenness with edge lengths
#'   \code{1 / weight}, and \code{"kmedoids"} selects medoids on the largest
#'   connected component using the same edge costs (farthest-first, then
#'   restricted PAM).
#' @param seed Random seed used when \code{method="random"} or
#'   \code{method="kmedoids"} (farthest-first start).
#' @param betweenness_cutoff Path-length cutoff for \code{igraph::betweenness}
#'   when \code{method="betweenness"}, measured on the sum of edge costs
#'   (\code{1/weight}). Default \code{-1} is exact betweenness. A positive value
#'   drops paths longer than that cost. Ignored for other methods.
#' @param kmedoids_m Per-cluster candidate cap (strength Top-M) for restricted
#'   PAM when \code{method="kmedoids"}. Ignored for other methods.
#' @param kmedoids_max_iter Maximum PAM passes when \code{method="kmedoids"}.
#'   Ignored for other methods.
#'
#' @return A character vector of landmark node IDs.
#' @examples
#' # Build a small graph
#' edges <- data.frame(
#'   node1 = c("A", "B", "C", "D", "E"),
#'   node2 = c("B", "C", "D", "E", "A"),
#'   weight = c(1.0, 2.0, 0.5, 1.5, 1.0)
#' )
#' adj <- gsemb_build_graph(edges, weight = "weight")
#' gsemb_select_landmarks(adj, k = 2, method = "degree")
#' gsemb_select_landmarks(adj, k = 2, method = "betweenness")
#'
#' @export
gsemb_select_landmarks <- function(adj,
                                   k = 128,
                                   method = c("degree", "random", "betweenness", "kmedoids"),
                                   seed = 1,
                                   betweenness_cutoff = -1,
                                   kmedoids_m = 50,
                                   kmedoids_max_iter = 10) {
  method <- match.arg(method)
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have rownames")
  k <- min(k, length(nodes))
  if (method %in% c("betweenness", "kmedoids")) {
    weights <- if (inherits(adj, "sparseMatrix")) adj@x else as.numeric(adj)
    if (any(!is.finite(weights) | weights < 0)) {
      stop("nonzero edge weights must be finite and positive")
    }
  }
  switch(
    method,
    degree = .select_landmarks_by_score(nodes, Matrix::rowSums(adj), k),
    random = {
      set.seed(seed)
      sample(nodes, k)
    },
    betweenness = .select_landmarks_betweenness(
      adj, nodes, k, betweenness_cutoff
    ),
    kmedoids = .select_landmarks_kmedoids(
      adj, k, seed, kmedoids_m, kmedoids_max_iter
    )
  )
}

.select_landmarks_by_score <- function(nodes, score, k) {
  nodes[order(score, decreasing = TRUE)][seq_len(k)]
}

.select_landmarks_betweenness <- function(adj, nodes, k, cutoff) {
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  costs <- 1 / pmax(igraph::E(g)$weight, .Machine$double.eps)
  bc <- igraph::betweenness(
    g,
    directed = FALSE,
    weights = costs,
    cutoff = cutoff
  )
  .select_landmarks_by_score(nodes, bc, k)
}

.select_landmarks_kmedoids <- function(adj, k, seed, m, max_iter) {
  g <- igraph::graph_from_adjacency_matrix(
    adj != 0,
    mode = "undirected",
    diag = FALSE
  )
  comp <- igraph::components(g)
  keep <- which(comp$membership == which.max(comp$csize))
  nodes <- igraph::V(g)$name[keep]
  n <- length(nodes)
  k <- min(as.integer(k), n)
  if (k >= n) {
    return(nodes)
  }

  # CSC of the induced subgraph for C++ SSSP / restricted PAM.
  sub <- Matrix::drop0(adj[nodes, nodes, drop = FALSE])
  idx0 <- kmedoids_hop(
    as.integer(sub@i),
    as.integer(sub@p),
    as.numeric(sub@x),
    as.integer(k),
    as.integer(seed),
    as.integer(m),
    as.integer(max_iter)
  )
  nodes[as.integer(idx0) + 1L]
}
