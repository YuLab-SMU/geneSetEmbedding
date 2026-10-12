gsemb_compute_diffusion_distributions <- function(adj,
                                                  seed_matrix,
                                                  alpha = 0.5,
                                                  tol = 1e-10,
                                                  max_iter = 200,
                                                  normalize = "col") {
  W <- gsemb_transition_matrix(adj, normalize = normalize)
  P <- gsemb_diffuse_seeds(W, seed_matrix, alpha = alpha, tol = tol, max_iter = max_iter)
  rownames(P) <- rownames(adj)
  P
}

#' Train Set2Gaussian embeddings with torch
#'
#' Learn a gene embedding (points) and gene set embeddings (diagonal Gaussians)
#' by matching diffusion distributions from landmarks and gene sets.
#'
#' @param adj Adjacency matrix with node IDs in rownames.
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param landmarks Optional landmark IDs. If `NULL`, selected with
#'   `landmark_method`.
#' @param k Number of landmarks when `landmarks` is `NULL`.
#' @param dim Embedding dimension.
#' @param alpha Restart probability for diffusion.
#' @param tol Convergence tolerance for diffusion.
#' @param max_iter Maximum diffusion iterations.
#' @param normalize Normalization for `gsemb_transition_matrix`.
#' @param tau Softmax temperature used in distribution matching.
#' @param lambda_set Weight for the set-level loss term.
#' @param epochs Number of optimization epochs.
#' @param lr Learning rate.
#' @param batch_size Batch size for sampling landmarks/sets per epoch.
#' @param seed Random seed.
#' @param landmark_method Landmark selection method when `landmarks` is `NULL`.
#' @param betweenness_cutoff,kmedoids_m,kmedoids_max_iter Passed to
#'   [gsemb_select_landmarks()] when `landmarks` is `NULL`.
#' @param device `"cpu"` or `"cuda"` (falls back to CPU if CUDA is unavailable).
#'
#' @return A list with `gene_embedding`, `set_mu`, `set_var`, `landmarks`, and `losses`.
#' @examples
#' \dontrun{
#' # Simulate a small adjacency matrix and gene sets
#' set.seed(42)
#' nodes <- paste0("GENE", 1:30)
#' edges <- data.frame(
#'   node1 = sample(nodes, 100, replace = TRUE),
#'   node2 = sample(nodes, 100, replace = TRUE)
#' )
#' adj <- gsemb_build_graph(edges)
#'
#' gene_sets <- list(
#'   SET1 = sample(nodes, 10),
#'   SET2 = sample(nodes, 8)
#' )
#'
#' # Fit Set2Gaussian model (requires torch)
#' fit <- gsemb_fit_set2gaussian_torch(
#'   adj,
#'   gene_sets,
#'   k = 10,
#'   dim = 16,
#'   alpha = 0.5,
#'   epochs = 5,
#'   lr = 5e-3,
#'   batch_size = 4,
#'   device = "cpu"
#' )
#'
#' # Inspect results
#' dim(fit$gene_embedding)
#' dim(fit$set_mu)
#' }
#' @export
gsemb_fit_set2gaussian_torch <- function(adj,
                                         gene_sets,
                                         landmarks = NULL,
                                         k = 128,
                                         dim = 64,
                                         alpha = 0.5,
                                         tol = 1e-10,
                                         max_iter = 200,
                                         normalize = "col",
                                         tau = 1.0,
                                         lambda_set = 1.0,
                                         epochs = 200,
                                         lr = 5e-3,
                                         batch_size = 8,
                                         seed = 1,
                                         landmark_method = c("degree", "random", "betweenness", "kmedoids"),
                                         betweenness_cutoff = -1,
                                         kmedoids_m = 50,
                                         kmedoids_max_iter = 10,
                                         device = c("cpu", "cuda")) {
  require_torch()
  device <- match.arg(device)
  landmark_method <- match.arg(landmark_method)
  gene_sets <- validate_gene_sets(gene_sets)
  if (!inherits(adj, "Matrix")) stop("adj must be a Matrix")
  nodes <- rownames(adj)
  if (is.null(nodes)) stop("adj must have rownames")
  n <- length(nodes)

  if (is.null(landmarks)) {
    landmarks <- gsemb_select_landmarks(
      adj,
      k = k,
      method = landmark_method,
      seed = seed,
      betweenness_cutoff = betweenness_cutoff,
      kmedoids_m = kmedoids_m,
      kmedoids_max_iter = kmedoids_max_iter
    )
  }
  lm_idx <- match(landmarks, nodes)
  if (any(is.na(lm_idx))) stop("some landmarks are not in graph nodes")

  S_lm <- matrix(0, n, length(lm_idx))
  S_lm[cbind(lm_idx, seq_along(lm_idx))] <- 1
  P_lm <- gsemb_compute_diffusion_distributions(adj, S_lm, alpha = alpha, tol = tol, max_iter = max_iter, normalize = normalize)
  P_lm <- t(P_lm)

  set_ids <- names(gene_sets)
  S_set <- matrix(0, n, length(set_ids))
  for (j in seq_along(set_ids)) {
    genes <- intersect(gene_sets[[set_ids[j]]], nodes)
    if (length(genes) == 0) next
    idx <- match(genes, nodes)
    S_set[idx, j] <- 1 / length(idx)
  }
  keep_sets <- colSums(S_set) > 0
  set_ids <- set_ids[keep_sets]
  S_set <- S_set[, keep_sets, drop = FALSE]
  P_set <- gsemb_compute_diffusion_distributions(adj, S_set, alpha = alpha, tol = tol, max_iter = max_iter, normalize = normalize)
  P_set <- t(P_set)

  # torch <- getNamespace("torch")
  # set.seed(seed)
  torch::torch_manual_seed(seed)
  dev <- if (device == "cuda" && torch::cuda_is_available()) torch::torch_device("cuda") else torch::torch_device("cpu")

  Z <- torch::torch_randn(c(n, dim), device = dev, dtype = torch::torch_float(), requires_grad = TRUE)
  mu <- torch::torch_randn(c(length(set_ids), dim), device = dev, dtype = torch::torch_float(), requires_grad = TRUE)
  logvar <- torch::torch_full(c(length(set_ids), dim), -2, device = dev, dtype = torch::torch_float(), requires_grad = TRUE)

  opt <- torch::optim_adam(list(Z, mu, logvar), lr = lr)

  target_lm <- torch::torch_tensor(P_lm, device = dev, dtype = torch::torch_float())
  target_set <- torch::torch_tensor(P_set, device = dev, dtype = torch::torch_float())

  z_norm2 <- function(x) {
    torch::torch_sum(x * x, dim = 2, keepdim = TRUE)
  }

  loss_landmark <- function(seed_idx, target_probs) {
    z_seed <- Z[seed_idx, , drop = FALSE]
    a2 <- z_norm2(z_seed)
    b2 <- z_norm2(Z)$transpose(1, 2)
    logits <- -(a2 + b2 - 2 * torch::torch_matmul(z_seed, Z$transpose(1, 2))) / tau
    log_q <- torch::nnf_log_softmax(logits, dim = 2)
    -(target_probs * log_q)$sum(dim = 2)$mean()
  }

  loss_set <- function(mu_batch, logvar_batch, target_probs) {
    var_batch <- torch::torch_exp(logvar_batch)
    z <- Z$unsqueeze(1)
    mu_b <- mu_batch$unsqueeze(2)
    var_b <- var_batch$unsqueeze(2)
    d2 <- ((z - mu_b)$pow(2) / var_b)$sum(dim = 3)
    cst <- logvar_batch$sum(dim = 2)
    logits <- -0.5 * (d2 + cst$unsqueeze(2))
    logits <- logits / tau
    log_q <- torch::nnf_log_softmax(logits, dim = 2)
    -(target_probs * log_q)$sum(dim = 2)$mean()
  }

  losses <- numeric(epochs)
  lm_all <- seq_len(length(lm_idx))
  set_all <- seq_len(length(set_ids))

  for (ep in seq_len(epochs)) {
    opt$zero_grad()

    lm_batch <- sample(lm_all, min(batch_size, length(lm_all)))
    lm_seed_idx <- lm_idx[lm_batch]
    l1 <- loss_landmark(lm_seed_idx, target_lm[lm_batch, , drop = FALSE])

    set_batch <- sample(set_all, min(batch_size, length(set_all)))
    l2 <- loss_set(mu[set_batch, , drop = FALSE], logvar[set_batch, , drop = FALSE], target_set[set_batch, , drop = FALSE])

    loss <- l1 + lambda_set * l2
    loss$backward()
    opt$step()

    losses[ep] <- as.numeric(loss$detach()$to(device = dev))
  }

  Z_cpu <- Z$detach()$to(device = dev)$to(dtype = torch::torch_double())
  mu_cpu <- mu$detach()$to(device = dev)$to(dtype = torch::torch_double())
  logvar_cpu <- logvar$detach()$to(device = dev)$to(dtype = torch::torch_double())

  gene_emb <- as.matrix(Z_cpu)
  rownames(gene_emb) <- nodes
  colnames(gene_emb) <- paste0("V", seq_len(dim))

  set_mu <- as.matrix(mu_cpu)
  rownames(set_mu) <- set_ids
  colnames(set_mu) <- paste0("V", seq_len(dim))

  set_var <- as.matrix(torch::torch_exp(logvar_cpu))
  rownames(set_var) <- set_ids
  colnames(set_var) <- paste0("V", seq_len(dim))

  list(
    gene_embedding = gene_emb,
    set_mu = set_mu,
    set_var = set_var,
    landmarks = landmarks,
    losses = losses,
    device = as.character(dev$type)
  )
}
