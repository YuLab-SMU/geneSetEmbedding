library(testthat)
library(geneSetEmbedding)

test_that("gsemb_fit_set2gaussian_torch runs when torch is available", {
    skip_if_not_installed("torch")

    edges <- data.frame(
        node1 = c("A", "A", "B", "C", "D"),
        node2 = c("B", "C", "C", "D", "E"),
        weight = c(1, 1, 1, 1, 1)
    )
    adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)
    gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("D", "E"))

    fit <- gsemb_fit_set2gaussian_torch(
        adj = adj,
        gene_sets = gene_sets,
        k = 2,
        dim = 4,
        epochs = 2,
        batch_size = 2,
        max_iter = 30,
        lr = 1e-2,
        device = "cpu"
    )

    expect_true(is.matrix(fit$gene_embedding))
    expect_true(is.matrix(fit$set_mu))
    expect_true(is.matrix(fit$set_var))
    expect_equal(nrow(fit$gene_embedding), 5)
    expect_equal(ncol(fit$gene_embedding), 4)
    expect_equal(nrow(fit$set_mu), 2)
    expect_equal(ncol(fit$set_mu), 4)
    expect_equal(length(fit$losses), 2)
})
