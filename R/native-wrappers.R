#' Diagonal 2-Wasserstein distance (fully vectorized via BLAS)
#'
#' @param mu n1 x d matrix of means
#' @param var n1 x d matrix of variances (must be > 0)
#' @param mu2 n2 x d matrix of means
#' @param var2 n2 x d matrix of variances (must be > 0)
#' @return n1 x n2 matrix of squared distances
#' @keywords internal
w2_distance <- function(mu, var, mu2, var2) {
  .Call(`_geneSetEmbedding_w2_distance`, mu, var, mu2, var2)
}

#' Symmetric KL divergence between diagonal Gaussians
#'
#' @param mu n1 x d matrix of means
#' @param var n1 x d matrix of variances
#' @param mu2 n2 x d matrix of means
#' @param var2 n2 x d matrix of variances
#' @return n1 x n2 matrix of symmetric KL divergences
#' @keywords internal
sym_kl_distance <- function(mu, var, mu2, var2) {
  .Call(`_geneSetEmbedding_sym_kl_distance`, mu, var, mu2, var2)
}

#' Parallel permutations for gene-set enrichment
#'
#' @param W n_genes x n_sets weight matrix
#' @param stats n_genes gene-level statistics
#' @param n_perm number of permutations
#' @param seed random seed
#' @return n_perm x n_sets matrix of null enrichment scores
#' @keywords internal
rcpp_enrichment_permutations <- function(W, stats, n_perm, seed) {
  .Call(`_geneSetEmbedding_rcpp_enrichment_permutations`, W, stats, n_perm, seed)
}

#' @keywords internal
#' @noRd
kmedoids_hop <- function(i, p, x, k, seed, m, max_iter) {
  .Call(
    `_geneSetEmbedding_kmedoids_hop`,
    i, p, x, k, seed, m, max_iter
  )
}
