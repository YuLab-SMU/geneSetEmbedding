library(testthat)
library(geneSetEmbedding)

test_that("Graph, diffusion, embedding, and concise sets run", {
    edges <- data.frame(
        node1 = c("A", "A", "B", "C", "D"),
        node2 = c("B", "C", "C", "D", "E"),
        weight = c(1, 1, 1, 1, 1)
    )
    adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)
    expect_true(inherits(adj, "dgCMatrix"))

    landmarks <- gsemb_select_landmarks(adj, k = 2, method = "degree")
    feat <- gsemb_compute_node_landmark_features(adj, landmarks = landmarks, alpha = 0.5, max_iter = 50)
    expect_equal(nrow(feat), 5)
    expect_equal(ncol(feat), 2)

    emb_res <- gsemb_fit_gene_embedding(feat, dim = 2, method = "svd")
    gene_emb <- emb_res$embedding
    expect_equal(dim(gene_emb), c(5, 2))

    gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("D", "E"))
    gauss <- gsemb_fit_set_gaussians_from_members(gene_emb, gene_sets)
    expect_equal(dim(gauss$mu), c(2, 2))
    expect_equal(dim(gauss$var), c(2, 2))

    sim <- gsemb_gene_cosine_similarity(gene_emb)
    expect_equal(dim(sim), c(5, 5))
    sim2 <- gsemb_row_cosine_similarity(gene_emb)
    expect_equal(dim(sim2), c(5, 5))

    d <- gsemb_set_gaussian_distance(gauss$mu, gauss$var, metric = "w2")
    expect_equal(dim(d), c(2, 2))
    d2 <- gsemb_diag_gaussian_w2(gauss$mu, gauss$var)
    expect_equal(dim(d2), c(2, 2))

    concise <- gsemb_make_concise_gene_sets(gene_emb, gauss$mu, gauss$var, gene_sets = gene_sets, top_n = 2, min_size = 1, restrict_to_members = TRUE)
    expect_equal(length(concise$S1), 2)
    concise0 <- gsemb_make_concise_gene_sets_from_gaussian(gene_emb, gauss$mu, gauss$var, gene_sets = gene_sets, top_n = 2, min_size = 1, restrict_to_members = TRUE)
    expect_equal(length(concise0$S1), 2)

    concise2 <- gsemb_make_concise_gene_sets(gene_emb, gauss$mu, gauss$var, gene_sets = gene_sets, select = "softmax_mass", mass = 0.8, temperature = 1, min_size = 1, restrict_to_members = TRUE)
    expect_true(length(concise2$S1) >= 1)
})

test_that("High-level API runs end-to-end on CPU", {
    edges <- data.frame(
        node1 = c("A", "A", "B", "C", "D"),
        node2 = c("B", "C", "C", "D", "E"),
        weight = c(1, 1, 1, 1, 1)
    )
    gene_sets <- list(S1 = c("A", "B", "C"), S2 = c("D", "E"))

    fit <- gsemb_fit(edges, gene_sets, weight = "weight", method = "svd", dim = 2, k = 2, max_iter = 30, epochs = 3)
    expect_true(inherits(fit, "gsemb_embedding"))
    expect_equal(dim(fit$gene_embedding), c(5, 2))

    gs <- gsemb_gene_similarity(fit, genes = c("A", "B"), other_genes = c("C", "D"))
    expect_equal(dim(gs), c(2, 2))

    ss <- gsemb_set_similarity(fit, metric = "w2")
    expect_equal(dim(ss), c(2, 2))

    stats <- c(A = 2, B = 1, C = 0.5, D = -1, E = -2)
    ea <- geneSetEmbedding:::gsemb_embedding_enrichment(stats, fit, gene_sets = gene_sets, nperm = 50, seed = 1, top_genes = 3)
    expect_true(is.data.frame(ea))
    expect_true(all(c("ID", "ES", "pvalue", "p.adjust", "core_enrichment") %in% names(ea)))
    expect_equal(nrow(ea), 2)

    concise <- gsemb_concise_gene_sets(fit, gene_sets = gene_sets, top_n = 2, min_size = 1, restrict_to_members = TRUE)
    expect_equal(length(concise$S1), 2)

    clusters <- list(C1 = c("A", "B", "C"), C2 = c("D", "E"))
    cs <- gsemb_cluster_similarity(fit, clusters, metric = "w2", min_size = 2)
    expect_equal(dim(cs), c(2, 2))
})
