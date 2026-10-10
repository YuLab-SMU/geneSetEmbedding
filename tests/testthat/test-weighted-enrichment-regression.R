library(testthat)

test_that("Weighted ORA and GSEA align genes and respect member restrictions", {
  E <- matrix(c(-1, 0, 1, 2), 4, 1,
              dimnames = list(c("A","B","C","D"), "V1"))
  mu <- matrix(0, 1, 1, dimnames = list("S1", "V1"))
  v <- matrix(1, 1, 1, dimnames = dimnames(mu))
  fit <- structure(list(gene_embedding = E, set_mu = mu, set_var = v),
                   class = "gsemb_embedding")
  y <- c(A = 0, B = 1, C = 0, D = 0)
  sets <- list(S1 = c("A","B"))
  ora <- gsemb_weighted_ora(y, fit, gene_sets = sets, nperm = 0,
                          score = "neg_mahalanobis", top_genes = 0)
  gsea <- gsemb_weighted_gsea(y[c("D","B","A","C")], fit,
                            gene_sets = sets, nperm = 0,
                            score = "neg_mahalanobis", top_genes = 0)
  expected <- 1 / (1 + 2 * exp(-1) + exp(-4))
  expect_equal(ora$ES, expected, tolerance = 1e-12)
  expect_equal(gsea$ES, expected, tolerance = 1e-12)
  restricted <- gsemb_weighted_ora(
    y, fit, gene_sets = sets, nperm = 0, score = "neg_mahalanobis",
    restrict_to_members = TRUE, top_genes = 0
  )
  expect_equal(restricted$ES, 1 / (1 + exp(-1)), tolerance = 1e-12)
  expect_error(gsemb_weighted_ora(c(A = 0, B = 0.5, C = 1, D = 0), fit),
               "0")
})

test_that("Monte Carlo p-values use the add-one correction", {
  E <- matrix(c(-1, 0, 1, 2), 4, 1,
              dimnames = list(c("A","B","C","D"), "V1"))
  mu <- matrix(0, 1, 1, dimnames = list("S1", "V1"))
  v <- matrix(1, 1, 1, dimnames = dimnames(mu))
  batch <- get(".gsemb_weighted_enrichment_sets",
               envir = environment(gsemb_weighted_gsea))
  out <- batch(
    gene_stats = c(A = 0, B = 1, C = 0, D = 0),
    gene_emb = E, genes = rownames(E), mu = mu, var = v,
    gene_sets = list(S1 = rownames(E)), restrict_to_members = FALSE,
    score = "neg_mahalanobis", temperature = 1,
    nperm = 4, alternative = "greater", seed = 1, eps = 1e-8,
    top_genes = 0, perm_mat = diag(4)
  )
  # Four possible positions of the binary label; only B reaches the observed score.
  expect_equal(out$pvalue, 2/5, tolerance = 1e-12)
  expect_equal(out$status, "ok")
})

test_that("Hub beta preserves identity at zero and reweights by node strength", {
  reweight <- get(".gsemb_apply_degree_beta",
                  envir = environment(gsemb_weighted_gsea))
  W <- matrix(c(0.5, 0.5), 2, 1, dimnames = list(c("A","B"), "S1"))
  expect_identical(reweight(W, c(A = 1, B = 4), beta = 0), W)
  expect_equal(unname(reweight(W, c(A = 1, B = 4), beta = 1)),
               matrix(c(0.8, 0.2), 2, 1), tolerance = 1e-12)
})
