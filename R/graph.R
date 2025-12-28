#' Build a PPI graph adjacency matrix
#'
#' Construct a sparse adjacency matrix from an edge list. Node names are stored
#' in the matrix dimnames and are used throughout the package.
#'
#' @param edges A data.frame containing at least two columns for endpoints.
#' @param node1,node2 Column names in \code{edges} for source/target.
#' @param weight Optional column name in \code{edges} providing edge weights.
#' @param directed Logical; whether to keep the graph directed.
#' @param nodes Optional character vector of node IDs to keep/order.
#'
#' @return A sparse adjacency matrix of class \code{dgCMatrix}.
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
#' Choose landmark nodes either by degree (highest first) or random sampling.
#'
#' @param adj Adjacency matrix with rownames as node IDs.
#' @param k Number of landmarks to select (capped at number of nodes).
#' @param method \code{"degree"} or \code{"random"}.
#' @param seed Random seed used when \code{method="random"}.
#'
#' @return A character vector of landmark node IDs.
#' @export
gsemb_select_landmarks <- function(adj, k = 128, method = c("degree", "random"), seed = 1) {
    method <- match.arg(method)
    nodes <- rownames(adj)
    if (is.null(nodes)) stop("adj must have rownames")
    n <- length(nodes)
    k <- min(k, n)
    if (method == "degree") {
        deg <- Matrix::rowSums(adj != 0)
        ord <- order(deg, decreasing = TRUE)
        nodes[ord][seq_len(k)]
    } else {
        set.seed(seed)
        sample(nodes, k)
    }
}
