#' Fit gene and gene set embeddings from PPI and gene sets
#'
#' This is the main entry point: build a graph from a PPI edge list (or accept
#' a pre-built adjacency matrix), compute diffusion features, learn a gene
#' embedding, and derive diagonal Gaussian embeddings for gene sets. Optionally,
#' train a Set2Gaussian model with torch.
#'
#' @param ppi Either a data.frame edge list or a Matrix adjacency matrix.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param node1,node2 Column names in `ppi` when `ppi` is a data.frame.
#' @param weight Optional weight column name in `ppi` when `ppi` is a data.frame.
#' @param directed Logical; whether the graph is directed.
#' @param nodes Optional node ID vector used to restrict/order the graph.
#' @param method Embedding method: `"svd"`, `"torch_autoencoder"`, or `"set2gaussian_torch"`.
#' @param dim Embedding dimension.
#' @param k Number of landmarks used for diffusion features.
#' @param alpha Restart probability for diffusion.
#' @param tol Convergence tolerance for diffusion.
#' @param max_iter Maximum diffusion iterations.
#' @param normalize Normalization for `gsemb_transition_matrix`.
#' @param epochs,lr,batch_size Training hyperparameters for torch methods.
#' @param seed Random seed.
#' @param device `"cpu"` or `"cuda"` (falls back to CPU if CUDA is unavailable).
#' @param ... Passed through to the selected training backend.
#'
#' @return A `gsemb_embedding` object with components:
#'   `gene_embedding`, `set_mu`, `set_var`, `adj`, and training metadata.
#' @examples
#' \dontrun{
#' # Simulate a small PPI edge list
#' set.seed(42)
#' nodes <- paste0("GENE", 1:50)
#' edges <- data.frame(
#'   node1 = sample(nodes, 200, replace = TRUE),
#'   node2 = sample(nodes, 200, replace = TRUE),
#'   weight = runif(200, 0.5, 1)
#' )
#' # Remove self-loops
#' edges <- edges[edges$node1 != edges$node2, ]
#'
#' # Simulate a few gene sets
#' gene_sets <- list(
#'   SET1 = sample(nodes, 10),
#'   SET2 = sample(nodes, 8),
#'   SET3 = sample(nodes, 12)
#' )
#'
#' # Fit embedding using SVD (fastest method)
#' fit <- gsemb_fit(
#'   ppi = edges,
#'   gene_sets = gene_sets,
#'   method = "svd",
#'   dim = 16,
#'   k = 20,
#'   alpha = 0.5,
#'   epochs = 5
#' )
#'
#' # Inspect the result
#' print(fit)
#'
#' # Compute gene-gene similarities
#' sim <- gsemb_gene_similarity(fit)
#' dim(sim)
#' }
#' @export
gsemb_fit <- function(ppi,
                      gene_sets,
                      node1 = "node1",
                      node2 = "node2",
                      weight = NULL,
                      directed = FALSE,
                      nodes = NULL,
                      method = c("svd", "torch_autoencoder", "set2gaussian_torch"),
                      dim = 64,
                      k = 128,
                      alpha = 0.5,
                      tol = 1e-10,
                      max_iter = 200,
                      normalize = "col",
                      epochs = 200,
                      lr = 5e-3,
                      batch_size = 8,
                      seed = 1,
                      device = c("cpu", "cuda"),
                      ...) {
  method <- match.arg(method)
  device <- match.arg(device)
  gene_sets <- validate_gene_sets(gene_sets)

  adj <- if (inherits(ppi, "Matrix")) {
    ppi
  } else {
    gsemb_build_graph(ppi, node1 = node1, node2 = node2, weight = weight, directed = directed, nodes = nodes)
  }
  if (is.null(rownames(adj))) stop("PPI graph must have node names")

  if (method == "set2gaussian_torch") {
    fit <- gsemb_fit_set2gaussian_torch(
      adj = adj,
      gene_sets = gene_sets,
      k = k,
      dim = dim,
      alpha = alpha,
      tol = tol,
      max_iter = max_iter,
      normalize = normalize,
      epochs = epochs,
      lr = lr,
      batch_size = batch_size,
      seed = seed,
      device = device,
      ...
    )
    res <- list(
      adj = adj,
      method = method,
      gene_embedding = fit$gene_embedding,
      set_mu = fit$set_mu,
      set_var = fit$set_var,
      landmarks = fit$landmarks,
      losses = fit$losses
    )
    class(res) <- "gsemb_embedding"
    return(res)
  }

  landmarks <- gsemb_select_landmarks(adj, k = k, method = "degree", seed = seed)
  node_features <- gsemb_compute_node_landmark_features(
    adj = adj,
    landmarks = landmarks,
    alpha = alpha,
    tol = tol,
    max_iter = max_iter,
    normalize = normalize,
    seed = seed
  )
  emb_fit <- gsemb_fit_gene_embedding(
    node_features = node_features,
    dim = dim,
    method = if (method == "torch_autoencoder") "torch_autoencoder" else "svd",
    epochs = epochs,
    lr = lr,
    batch_size = max(64, 4 * dim),
    seed = seed,
    device = device
  )
  gene_embedding <- emb_fit$embedding
  gauss <- gsemb_fit_set_gaussians_from_members(gene_embedding, gene_sets)

  res <- list(
    adj = adj,
    method = method,
    gene_embedding = gene_embedding,
    set_mu = gauss$mu,
    set_var = gauss$var,
    landmarks = landmarks,
    losses = NULL
  )
  class(res) <- "gsemb_embedding"
  res
}

#' Compute gene-gene cosine similarities
#'
#' @param x A `gsemb_embedding` object or a numeric gene embedding matrix.
#' @param genes Optional character vector of genes to include as rows.
#' @param other_genes Optional character vector of genes to include as columns.
#' @param eps Small constant for numerical stability.
#'
#' @return A numeric similarity matrix.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' sim_all <- gsemb_gene_similarity(fit)
#'
#' # Compute similarities for a subset of genes
#' sim_sub <- gsemb_gene_similarity(fit, genes = c("GENE1", "GENE2", "GENE3"))
#' }
#' @export
gsemb_gene_similarity <- function(x,
                                  genes = NULL,
                                  other_genes = NULL,
                                  eps = 1e-12) {
  gene_emb <- if (inherits(x, "gsemb_embedding")) x$gene_embedding else x
  if (is.null(genes) && is.null(other_genes)) {
    return(gsemb_gene_cosine_similarity(gene_emb, eps = eps))
  }
  all_genes <- rownames(gene_emb)
  if (is.null(all_genes)) stop("gene embedding must have rownames")
  if (is.null(genes)) genes <- all_genes
  genes <- intersect(as.character(genes), all_genes)
  if (length(genes) == 0) stop("no genes found in embedding")
  if (is.null(other_genes)) {
    other_genes <- genes
  } else {
    other_genes <- intersect(as.character(other_genes), all_genes)
    if (length(other_genes) == 0) stop("no other_genes found in embedding")
  }
  gsemb_gene_cosine_similarity(gene_emb[genes, , drop = FALSE], gene_emb[other_genes, , drop = FALSE], eps = eps)
}

#' Compute set-set distances between Gaussian embeddings
#'
#' @param x A `gsemb_embedding` object.
#' @param sets Optional character vector of sets to include as rows.
#' @param other_sets Optional character vector of sets to include as columns.
#' @param metric Distance metric: `"w2"` or `"sym_kl"`.
#' @param eps Lower bound for variances.
#'
#' @return A numeric distance matrix.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Compute Wasserstein-2 distances between all gene sets
#' dist_w2 <- gsemb_set_similarity(fit, metric = "w2")
#'
#' # Compute symmetric KL divergences for a subset of sets
#' dist_kl <- gsemb_set_similarity(
#'   fit,
#'   sets = c("SET1", "SET2"),
#'   metric = "sym_kl"
#' )
#' }
#' @export
gsemb_set_similarity <- function(x,
                                 sets = NULL,
                                 other_sets = NULL,
                                 metric = c("w2", "sym_kl"),
                                 eps = 1e-8) {
  metric <- match.arg(metric)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  mu <- x$set_mu
  var <- x$set_var
  if (is.null(rownames(mu))) stop("set embedding must have rownames")

  all_sets <- rownames(mu)
  if (is.null(sets) && is.null(other_sets)) {
    return(gsemb_set_gaussian_distance(mu, var, metric = metric, eps = eps))
  }
  if (is.null(sets)) sets <- all_sets
  sets <- intersect(as.character(sets), all_sets)
  if (length(sets) == 0) stop("no sets found in embedding")
  if (is.null(other_sets)) {
    other_sets <- sets
  } else {
    other_sets <- intersect(as.character(other_sets), all_sets)
    if (length(other_sets) == 0) stop("no other_sets found in embedding")
  }
  gsemb_set_gaussian_distance(mu[sets, , drop = FALSE], var[sets, , drop = FALSE], mu[other_sets, , drop = FALSE], var[other_sets, , drop = FALSE], metric = metric, eps = eps)
}

#' Compute distances between two collections of gene clusters
#'
#' Gene clusters are represented as diagonal Gaussians fitted from member genes
#' in the learned gene embedding.
#'
#' @param x A `gsemb_embedding` object.
#' @param clusters Named list of character vectors (gene IDs).
#' @param other_clusters Optional second named list of clusters.
#' @param metric Distance metric: `"w2"` or `"sym_kl"`.
#' @param eps Lower bound for variances.
#' @param min_size Minimum number of genes per cluster after intersecting with the embedding.
#'
#' @return A numeric distance matrix (clusters x other_clusters).
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define some custom clusters
#' clusters <- list(
#'   CLUST1 = c("GENE1", "GENE2", "GENE3"),
#'   CLUST2 = c("GENE4", "GENE5", "GENE6")
#' )
#' # Compute distances between these clusters
#' dist_mat <- gsemb_cluster_similarity(fit, clusters = clusters, metric = "w2")
#'
#' # Compute distances between two different cluster lists
#' other_clusters <- list(
#'   OTHER1 = c("GENE10", "GENE11"),
#'   OTHER2 = c("GENE12", "GENE13", "GENE14")
#' )
#' dist_mat2 <- gsemb_cluster_similarity(
#'   fit,
#'   clusters = clusters,
#'   other_clusters = other_clusters,
#'   metric = "sym_kl"
#' )
#' }
#' @export
gsemb_cluster_similarity <- function(x,
                                     clusters,
                                     other_clusters = NULL,
                                     metric = c("w2", "sym_kl"),
                                     eps = 1e-8,
                                     min_size = 2) {
  metric <- match.arg(metric)
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  if (!is.list(clusters) || is.null(names(clusters))) stop("clusters must be a named list")
  if (is.null(other_clusters)) other_clusters <- clusters
  if (!is.list(other_clusters) || is.null(names(other_clusters))) stop("other_clusters must be a named list")

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) stop("gene embedding must have rownames")

  to_gauss <- function(gs) {
    gs <- lapply(gs, function(v) intersect(as.character(v), rownames(gene_emb)))
    gs <- gs[vapply(gs, length, integer(1)) >= min_size]
    if (length(gs) == 0) stop("no clusters have enough genes in embedding")
    g <- gsemb_fit_set_gaussians_from_members(gene_emb, gs, eps = eps)
    list(mu = g$mu, var = g$var)
  }

  a <- to_gauss(clusters)
  b <- to_gauss(other_clusters)
  gsemb_set_gaussian_distance(a$mu, a$var, b$mu, b$var, metric = metric, eps = eps)
}

#' Create concise gene sets from a fitted embedding
#'
#' Convenience wrapper around `gsemb_make_concise_gene_sets` using the embedding
#' stored in a `gsemb_embedding` object.
#'
#' @param x A `gsemb_embedding` object.
#' @param gene_sets Named list of gene sets (used when restricting candidates).
#' @param ... Passed to `gsemb_make_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define original gene sets
#' orig_sets <- list(
#'   PATH1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
#'   PATH2 = c("GENE3", "GENE4", "GENE6", "GENE7")
#' )
#' # Create concise versions
#' concise <- gsemb_concise_gene_sets(fit, gene_sets = orig_sets, top_n = 3)
#'
#' # View concise sets
#' str(concise)
#' }
#' @export
gsemb_concise_gene_sets <- function(x,
                                    gene_sets,
                                    ...) {
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  gene_sets <- validate_gene_sets(gene_sets)
  gsemb_make_concise_gene_sets(
    gene_embedding = x$gene_embedding,
    set_mu = x$set_mu,
    set_var = x$set_var,
    gene_sets = gene_sets,
    ...
  )
}

#' One-click training from PPI edges and gene sets
#'
#' @param ppi_edges PPI edge list as a data.frame.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param ... Passed to `gsemb_fit`.
#'
#' @return A `gsemb_embedding` object.
#' @examples
#' \dontrun{
#' # Simulate a small PPI edge list
#' set.seed(42)
#' nodes <- paste0("GENE", 1:50)
#' edges <- data.frame(
#'   node1 = sample(nodes, 200, replace = TRUE),
#'   node2 = sample(nodes, 200, replace = TRUE),
#'   weight = runif(200, 0.5, 1)
#' )
#' # Remove self-loops
#' edges <- edges[edges$node1 != edges$node2, ]
#'
#' # Simulate a few gene sets
#' gene_sets <- list(
#'   SET1 = sample(nodes, 10),
#'   SET2 = sample(nodes, 8),
#'   SET3 = sample(nodes, 12)
#' )
#'
#' # One-click training
#' fit <- gsemb_train_embedding_from_ppi_and_genesets(
#'   ppi_edges = edges,
#'   gene_sets = gene_sets,
#'   method = "svd",
#'   dim = 16,
#'   k = 20,
#'   alpha = 0.5
#' )
#'
#' # Inspect the result
#' print(fit)
#' }
#' @export
gsemb_train_embedding_from_ppi_and_genesets <- function(ppi_edges,
                                                        gene_sets,
                                                        ...) {
  gsemb_fit(ppi = ppi_edges, gene_sets = gene_sets, ...)
}

#' Compute all default similarity matrices from a fitted embedding
#'
#' @param x A `gsemb_embedding` object.
#' @param eps_gene Small constant for cosine similarity stability.
#' @param eps_set Lower bound for Gaussian variances.
#'
#' @return A list with `gene_gene`, `set_set_w2`, and `set_set_kl`.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Compute all similarity matrices
#' all_sims <- gsemb_calculate_all_similarities(fit)
#'
#' # Access individual matrices
#' gene_gene_sim <- all_sims$gene_gene
#' set_set_w2 <- all_sims$set_set_w2
#' set_set_kl <- all_sims$set_set_kl
#' }
#' @export
gsemb_calculate_all_similarities <- function(x,
                                             eps_gene = 1e-12,
                                             eps_set = 1e-8) {
  if (!inherits(x, "gsemb_embedding")) stop("x must be a gsemb_embedding object")
  list(
    gene_gene = gsemb_gene_similarity(x, eps = eps_gene),
    set_set_w2 = gsemb_set_similarity(x, metric = "w2", eps = eps_set),
    set_set_kl = gsemb_set_similarity(x, metric = "sym_kl", eps = eps_set)
  )
}

#' Create concise gene sets (high-level wrapper)
#'
#' @param x A `gsemb_embedding` object.
#' @param gene_sets Named list of gene sets.
#' @param ... Passed to `gsemb_concise_gene_sets`.
#'
#' @return A named list of concise gene sets.
#' @examples
#' \dontrun{
#' # Assuming you have a fitted gsemb_embedding object named 'fit'
#' # (see example in gsemb_fit)
#' # Define original gene sets
#' orig_sets <- list(
#'   PATH1 = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
#'   PATH2 = c("GENE3", "GENE4", "GENE6", "GENE7")
#' )
#' # Get concise versions using the wrapper
#' concise <- gsemb_get_concise_gene_sets(fit, gene_sets = orig_sets, top_n = 3)
#'
#' # View concise sets
#' str(concise)
#' }
#' @export
gsemb_get_concise_gene_sets <- function(x,
                                        gene_sets,
                                        ...) {
  gsemb_concise_gene_sets(x, gene_sets = gene_sets, ...)
}

# Soft-membership enrichment (Weighted-ORA / Weighted-GSEA, hub correction)
# lives in R/weighted_enrichment.R:
#   gsemb_embedding_enrichment(), gsemb_weighted_gsea(), gsemb_weighted_ora(),
#   gsemb_sweep_hub_beta().
