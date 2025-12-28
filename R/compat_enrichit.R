#' Coerce an embedding object to a numeric matrix
#'
#' @param x A numeric matrix, Matrix, or data.frame representing an embedding.
#' @param id_col Optional ID column name when `x` is a data.frame.
#'
#' @return A numeric matrix.
#' @export
gsemb_as_embedding_matrix <- function(x, id_col = NULL) {
    as_numeric_matrix(x, id_col = id_col)
}

#' Row-wise cosine similarity
#'
#' Compatibility wrapper for cosine similarity.
#'
#' @param x Numeric embedding matrix.
#' @param y Optional second embedding matrix.
#' @param eps Small constant for numerical stability.
#'
#' @return A cosine similarity matrix.
#' @export
gsemb_row_cosine_similarity <- function(x, y = NULL, eps = 1e-12) {
    gsemb_gene_cosine_similarity(gene_embedding = x, other = y, eps = eps)
}

#' Diagonal Gaussian 2-Wasserstein distance
#'
#' Compatibility wrapper for diagonal Gaussian distance with `metric="w2"`.
#'
#' @param mu Mean matrix for sets.
#' @param var Variance matrix for sets.
#' @param mu2 Optional mean matrix for other sets.
#' @param var2 Optional variance matrix for other sets.
#' @param eps Lower bound for variances.
#'
#' @return A distance matrix.
#' @export
gsemb_diag_gaussian_w2 <- function(mu, var, mu2 = NULL, var2 = NULL, eps = 1e-8) {
    gsemb_set_gaussian_distance(set_mu = mu, set_var = var, other_mu = mu2, other_var = var2, metric = "w2", eps = eps)
}

#' Diagonal Gaussian symmetric KL distance
#'
#' Compatibility wrapper for diagonal Gaussian distance with `metric="sym_kl"`.
#'
#' @param mu Mean matrix for sets.
#' @param var Variance matrix for sets.
#' @param mu2 Optional mean matrix for other sets.
#' @param var2 Optional variance matrix for other sets.
#' @param eps Lower bound for variances.
#'
#' @return A distance matrix.
#' @export
gsemb_diag_gaussian_sym_kl <- function(mu, var, mu2 = NULL, var2 = NULL, eps = 1e-8) {
    gsemb_set_gaussian_distance(set_mu = mu, set_var = var, other_mu = mu2, other_var = var2, metric = "sym_kl", eps = eps)
}

#' Gene-to-set score under diagonal Gaussian
#'
#' Compatibility wrapper for `gsemb_gene_to_set_score`.
#'
#' @param gene_emb Numeric gene embedding matrix.
#' @param set_mu Mean matrix for sets.
#' @param set_var Variance matrix for sets.
#' @param score Scoring function.
#' @param eps Lower bound for variances.
#'
#' @return A numeric matrix of scores.
#' @export
gsemb_gene_to_diag_gaussian_score <- function(gene_emb,
                                            set_mu,
                                            set_var,
                                            score = c("loglik", "neg_mahalanobis"),
                                            eps = 1e-8) {
    gsemb_gene_to_set_score(gene_embedding = gene_emb, set_mu = set_mu, set_var = set_var, score = score, eps = eps)
}

#' Make concise gene sets from diagonal Gaussian embeddings
#'
#' Compatibility wrapper for `gsemb_make_concise_gene_sets`.
#'
#' @param gene_emb Numeric gene embedding matrix.
#' @param set_mu Mean matrix for sets.
#' @param set_var Variance matrix for sets.
#' @param gene_sets Optional named list of gene sets used to restrict candidates.
#' @param top_n Number of genes per set when using `select="top_n"`.
#' @param restrict_to_members Restrict candidate genes to members of the input set.
#' @param min_size Minimum size of each output set.
#' @param max_size Maximum size of each output set.
#' @param score Scoring function.
#' @param eps Lower bound for variances.
#' @param ... Passed to `gsemb_make_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
#' @export
gsemb_make_concise_gene_sets_from_gaussian <- function(gene_emb,
                                                     set_mu,
                                                     set_var,
                                                     gene_sets = NULL,
                                                     top_n = 50,
                                                     restrict_to_members = TRUE,
                                                     min_size = 5,
                                                     max_size = 500,
                                                     score = c("loglik", "neg_mahalanobis"),
                                                     eps = 1e-8,
                                                     ...) {
    gsemb_make_concise_gene_sets(
        gene_embedding = gene_emb,
        set_mu = set_mu,
        set_var = set_var,
        gene_sets = gene_sets,
        top_n = top_n,
        restrict_to_members = restrict_to_members,
        min_size = min_size,
        max_size = max_size,
        score = score,
        eps = eps,
        ...
    )
}
