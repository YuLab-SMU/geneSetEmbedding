#' Random walk with restart (RWR)
#'
#' Compute a stationary diffusion distribution from a seed vector on a graph
#' transition matrix.
#'
#' @param W A transition matrix (typically from \code{gsemb_transition_matrix}).
#' @param seed Numeric vector of length \code{nrow(W)} giving restart mass.
#' @param alpha Restart probability.
#' @param tol Convergence tolerance on L1 change.
#' @param max_iter Maximum number of iterations.
#'
#' @return A numeric vector of length \code{nrow(W)} that sums to 1.
#' @export
gsemb_rwr <- function(W, seed, alpha = 0.5, tol = 1e-10, max_iter = 200) {
    if (!inherits(W, "Matrix")) stop("W must be a Matrix")
    seed <- as.numeric(seed)
    if (length(seed) != nrow(W)) stop("seed length must match nrow(W)")
    ssum <- sum(seed)
    if (ssum <= 0) stop("seed must have positive mass")
    seed <- seed / ssum

    p <- seed
    for (i in seq_len(max_iter)) {
        p_new <- as.numeric((1 - alpha) * (W %*% p) + alpha * seed)
        if (sum(abs(p_new - p)) < tol) {
            p <- p_new
            break
        }
        p <- p_new
    }
    p <- p / sum(p)
    p
}

#' Diffuse multiple seed vectors with RWR
#'
#' Vectorized RWR for a seed matrix where each column is a seed distribution.
#'
#' @param W A transition matrix (typically from \code{gsemb_transition_matrix}).
#' @param S Numeric seed matrix with \code{nrow(S) == nrow(W)}.
#' @param alpha Restart probability.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#'
#' @return A numeric matrix with the same dimension as \code{S}, each column summing to 1.
#' @export
gsemb_diffuse_seeds <- function(W, S, alpha = 0.5, tol = 1e-10, max_iter = 200) {
    if (!inherits(W, "Matrix")) stop("W must be a Matrix")
    if (!is.matrix(S)) stop("S must be a numeric matrix")
    if (nrow(S) != nrow(W)) stop("nrow(S) must match nrow(W)")
    col_sums <- colSums(S)
    if (any(col_sums <= 0)) stop("each seed column must have positive mass")
    S <- sweep(S, 2, col_sums, "/")

    P <- S
    for (i in seq_len(max_iter)) {
        P_new <- (1 - alpha) * (W %*% P) + alpha * S
        P_new <- as.matrix(P_new)
        if (sum(abs(P_new - P)) < tol * ncol(P)) {
            P <- P_new
            break
        }
        P <- P_new
    }
    P <- as.matrix(P)
    P <- sweep(P, 2, base::colSums(P), "/")
    P
}

#' Compute node diffusion features to landmark nodes
#'
#' For each landmark node, run diffusion with a one-hot seed and collect the
#' stationary distributions as node features.
#'
#' @param adj Adjacency matrix with node IDs in rownames.
#' @param landmarks Optional landmark node IDs; if \code{NULL} they are selected.
#' @param k Number of landmarks when \code{landmarks} is \code{NULL}.
#' @param alpha Restart probability.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param normalize Normalization for \code{gsemb_transition_matrix}.
#' @param landmark_method Landmark selection method when \code{landmarks} is \code{NULL}.
#' @param seed Random seed for landmark selection.
#'
#' @return A numeric matrix of size \code{n_nodes x n_landmarks}.
#' @export
gsemb_compute_node_landmark_features <- function(adj,
                                               landmarks = NULL,
                                               k = 128,
                                               alpha = 0.5,
                                               tol = 1e-10,
                                               max_iter = 200,
                                               normalize = "col",
                                               landmark_method = "degree",
                                               seed = 1) {
    if (is.null(landmarks)) {
        landmarks <- gsemb_select_landmarks(adj, k = k, method = landmark_method, seed = seed)
    }
    nodes <- rownames(adj)
    if (is.null(nodes)) stop("adj must have rownames")
    idx <- match(landmarks, nodes)
    if (any(is.na(idx))) stop("some landmarks are not in graph nodes")

    W <- gsemb_transition_matrix(adj, normalize = normalize)
    S <- matrix(0, nrow(adj), length(idx))
    S[cbind(idx, seq_along(idx))] <- 1
    P <- gsemb_diffuse_seeds(W, S, alpha = alpha, tol = tol, max_iter = max_iter)
    rownames(P) <- nodes
    colnames(P) <- landmarks
    P
}

#' Compute gene set diffusion features to landmark nodes
#'
#' For each gene set, diffuse from its member genes and extract the diffusion
#' mass on the provided landmark nodes.
#'
#' @param adj Adjacency matrix with node IDs in rownames.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param landmarks Character vector of landmark node IDs.
#' @param alpha Restart probability.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param normalize Normalization for \code{gsemb_transition_matrix}.
#' @param batch_size Number of gene sets to process per batch.
#'
#' @return A numeric matrix of size \code{n_sets x n_landmarks}.
#' @export
gsemb_compute_set_landmark_features <- function(adj,
                                              gene_sets,
                                              landmarks,
                                              alpha = 0.5,
                                              tol = 1e-10,
                                              max_iter = 200,
                                              normalize = "col",
                                              batch_size = 32) {
    gene_sets <- validate_gene_sets(gene_sets)
    nodes <- rownames(adj)
    if (is.null(nodes)) stop("adj must have rownames")
    W <- gsemb_transition_matrix(adj, normalize = normalize)
    lm_idx <- match(landmarks, nodes)
    if (any(is.na(lm_idx))) stop("some landmarks are not in graph nodes")

    set_ids <- names(gene_sets)
    out <- matrix(NA_real_, length(set_ids), length(landmarks))
    rownames(out) <- set_ids
    colnames(out) <- landmarks

    i <- 1
    while (i <= length(set_ids)) {
        j <- min(i + batch_size - 1, length(set_ids))
        batch_ids <- set_ids[i:j]
        S <- matrix(0, nrow(adj), length(batch_ids))
        for (b in seq_along(batch_ids)) {
            genes <- intersect(gene_sets[[batch_ids[b]]], nodes)
            if (length(genes) == 0) next
            idx <- match(genes, nodes)
            S[idx, b] <- 1 / length(idx)
        }
        keep <- colSums(S) > 0
        if (any(keep)) {
            P <- gsemb_diffuse_seeds(W, S[, keep, drop = FALSE], alpha = alpha, tol = tol, max_iter = max_iter)
            out[batch_ids[keep], ] <- t(P[lm_idx, , drop = FALSE])
        }
        i <- j + 1
    }
    out
}
