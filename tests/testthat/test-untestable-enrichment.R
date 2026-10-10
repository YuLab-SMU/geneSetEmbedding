library(testthat)
library(geneSetEmbedding)

test_that("Untestable sets get NA and are excluded from BH", {
  edges <- data.frame(
    node1 = c("A", "A", "B", "C", "D"),
    node2 = c("B", "C", "C", "D", "E"),
    weight = c(1, 1, 1, 1, 1)
  )
  gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("D", "E"))
  fit <- gsemb_fit(
    edges, gene_sets,
    weight = "weight", method = "svd", dim = 2, k = 2, max_iter = 30, epochs = 3
  )

  # Inject an empty-coverage pathway with NA Gaussian params (as fit would leave them).
  d <- ncol(fit$set_mu)
  fit$set_mu <- rbind(
    fit$set_mu,
    S_empty = rep(NA_real_, d)
  )
  fit$set_var <- rbind(
    fit$set_var,
    S_empty = rep(NA_real_, d)
  )
  gene_sets_ext <- c(gene_sets, list(S_empty = c("NOT_IN_GRAPH_1", "NOT_IN_GRAPH_2")))

  stats <- c(A = 2, B = 1, C = 0.5, D = -1, E = -2)
  res <- gsemb_weighted_gsea(
    stats, fit,
    gene_sets = gene_sets_ext,
    nperm = 50, seed = 1, top_genes = 3, alternative = "greater"
  )

  expect_equal(nrow(res), 3)
  expect_true("status" %in% names(res))

  empty_row <- res[res$ID == "S_empty", , drop = FALSE]
  expect_equal(nrow(empty_row), 1L)
  expect_true(empty_row$status %in% c("untestable_params", "untestable_coverage", "untestable_weights"))
  expect_true(is.na(empty_row$ES))
  expect_true(is.na(empty_row$z))
  expect_true(is.na(empty_row$pvalue))
  expect_true(is.na(empty_row$p.adjust))
  expect_true(is.na(empty_row$core_enrichment))

  ok <- res[res$status == "ok", , drop = FALSE]
  expect_true(nrow(ok) >= 1L)
  expect_true(all(is.finite(ok$pvalue)))
  expect_true(all(is.finite(ok$p.adjust)))

  # BH on ok subset only: matches p.adjust(ok$pvalue) (order-preserving).
  expect_equal(
    unname(ok$p.adjust),
    unname(stats::p.adjust(ok$pvalue, method = "BH")),
    tolerance = 1e-12
  )

  # Contaminating BH with a fake P=1 for the empty row would change ok adjustments
  # whenever there is more than one ok p-value (m differs).
  fake <- c(ok$pvalue, 1)
  fake_adj_ok <- stats::p.adjust(fake, method = "BH")[seq_len(nrow(ok))]
  if (nrow(ok) >= 2L) {
    expect_false(isTRUE(all.equal(unname(ok$p.adjust), unname(fake_adj_ok))))
  } else {
    # With a single ok test, BH(p) == p whether or not an extra 1 is appended;
    # still require the empty row itself is excluded (already checked NA above).
    expect_equal(unname(ok$p.adjust), unname(ok$pvalue), tolerance = 1e-12)
  }
})
