#' Compute cosine similarity between gene embeddings
#'
#' @param gene_embedding Numeric matrix with genes in rows and embedding dimensions in columns.
#'   Row names must be gene identifiers.
#' @param other Optional numeric matrix with the same number of columns as `gene_embedding`.
#'   If `NULL`, similarities are computed within `gene_embedding`.
#' @param eps Numeric scalar added to norms for numerical stability.
#'
#' @return A numeric matrix of cosine similarities. Row/column names are inherited from
#'   the row names of the input embedding matrices.
#'
#' @export
gsemb_gene_cosine_similarity <- function(gene_embedding, other = NULL, eps = 1e-12) {
    X <- as_numeric_matrix(gene_embedding)
    if (is.null(rownames(X))) stop("gene_embedding must have rownames")
    if (is.null(other)) {
        Y <- X
    } else {
        Y <- as_numeric_matrix(other)
        if (is.null(rownames(Y))) stop("other must have rownames")
    }
    if (ncol(X) != ncol(Y)) stop("embeddings must have the same number of columns")
    x_norm <- sqrt(rowSums(X * X))
    y_norm <- sqrt(rowSums(Y * Y))
    x_norm <- pmax(x_norm, eps)
    y_norm <- pmax(y_norm, eps)
    sim <- (X / x_norm) %*% t(Y / y_norm)
    rownames(sim) <- rownames(X)
    colnames(sim) <- rownames(Y)
    sim
}

#' Compute pairwise distances between diagonal Gaussian set embeddings
#'
#' @param set_mu Numeric matrix of Gaussian means (gene sets in rows).
#' @param set_var Numeric matrix of diagonal variances (same shape/rownames as `set_mu`).
#' @param other_mu Optional matrix of means for the second collection.
#' @param other_var Optional matrix of variances for the second collection.
#' @param metric Distance metric: `"w2"` for diagonal 2-Wasserstein distance, or `"sym_kl"`
#'   for symmetric KL divergence.
#' @param eps Numeric scalar lower bound for variances.
#'
#' @return A numeric matrix of pairwise distances.
#'
#' @export
gsemb_set_gaussian_distance <- function(set_mu,
                                      set_var,
                                      other_mu = NULL,
                                      other_var = NULL,
                                      metric = c("w2", "sym_kl"),
                                      eps = 1e-8) {
    metric <- match.arg(metric)
    mu <- as_numeric_matrix(set_mu)
    var <- as_numeric_matrix(set_var)
    if (is.null(rownames(mu)) || is.null(rownames(var))) stop("set_mu/set_var must have rownames")
    if (!all(dim(mu) == dim(var))) stop("set_mu and set_var must have identical dimensions")
    if (!all(rownames(mu) == rownames(var))) stop("set_mu and set_var must have identical rownames")

    if (is.null(other_mu)) other_mu <- mu
    if (is.null(other_var)) other_var <- var
    mu2 <- as_numeric_matrix(other_mu)
    var2 <- as_numeric_matrix(other_var)
    if (!all(dim(mu2) == dim(var2))) stop("other_mu and other_var must have identical dimensions")
    if (ncol(mu) != ncol(mu2)) stop("mu and other_mu must have same number of columns")

    var <- pmax(var, eps)
    var2 <- pmax(var2, eps)

    out <- matrix(0, nrow(mu), nrow(mu2))
    for (i in seq_len(nrow(mu))) {
        for (j in seq_len(nrow(mu2))) {
            dmu <- mu[i, ] - mu2[j, ]
            if (metric == "w2") {
                out[i, j] <- sum(dmu * dmu) + sum((sqrt(var[i, ]) - sqrt(var2[j, ]))^2)
            } else {
                kl01 <- 0.5 * (sum(log(var2[j, ] / var[i, ])) - ncol(mu) + sum(var[i, ] / var2[j, ]) + sum((dmu * dmu) / var2[j, ]))
                kl10 <- 0.5 * (sum(log(var[i, ] / var2[j, ])) - ncol(mu) + sum(var2[j, ] / var[i, ]) + sum((dmu * dmu) / var[i, ]))
                out[i, j] <- 0.5 * (kl01 + kl10)
            }
        }
    }
    rownames(out) <- rownames(mu)
    colnames(out) <- rownames(mu2)
    out
}
