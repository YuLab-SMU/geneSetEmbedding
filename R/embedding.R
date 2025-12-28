#' Fit a gene embedding from diffusion features
#'
#' Reduce node diffusion feature vectors into a lower-dimensional embedding using
#' either SVD (CPU) or a simple torch autoencoder.
#'
#' @param node_features Numeric matrix with nodes in rows and features in columns.
#'   Row names must be node IDs.
#' @param dim Embedding dimension.
#' @param method `"svd"` or `"torch_autoencoder"`.
#' @param epochs Number of epochs for the torch autoencoder.
#' @param lr Learning rate for the torch autoencoder.
#' @param batch_size Batch size for the torch autoencoder.
#' @param seed Random seed.
#' @param device `"cpu"` or `"cuda"` (falls back to CPU if CUDA is unavailable).
#'
#' @return A list with components `embedding` (numeric matrix) and `method`.
#'   When `method="torch_autoencoder"`, the fitted `model` is also returned.
#' @export
gsemb_fit_gene_embedding <- function(node_features,
                                   dim = 64,
                                   method = c("torch_autoencoder", "svd"),
                                   epochs = 50,
                                   lr = 1e-3,
                                   batch_size = 256,
                                   seed = 1,
                                   device = c("cpu", "cuda")) {
    method <- match.arg(method)
    device <- match.arg(device)

    X <- as_numeric_matrix(node_features)
    if (is.null(rownames(X))) stop("node_features must have rownames")
    if (dim <= 0 || dim > ncol(X)) stop("dim must be between 1 and ncol(node_features)")

    if (method == "svd") {
        Xc <- scale(X, center = TRUE, scale = FALSE)
        s <- svd(Xc, nu = dim, nv = 0)
        emb <- s$u %*% diag(s$d[seq_len(dim)], nrow = dim, ncol = dim)
        rownames(emb) <- rownames(X)
        colnames(emb) <- paste0("V", seq_len(dim))
        return(list(embedding = emb, model = NULL, method = method))
    }

    require_torch()
    torch <- getNamespace("torch")
    set.seed(seed)
    torch_manual_seed <- get("torch_manual_seed", envir = torch)
    torch_manual_seed(seed)

    dev <- if (device == "cuda" && torch::cuda_is_available()) torch::device("cuda") else torch::device("cpu")

    X_tensor <- torch::torch_tensor(X, dtype = torch::torch_float())$to(device = dev)

    encoder <- torch::nn_linear(ncol(X), dim)
    decoder <- torch::nn_linear(dim, ncol(X))
    model <- torch::nn_module(
        initialize = function() {
            self$enc <- encoder
            self$dec <- decoder
        },
        forward = function(x) {
            z <- self$enc(x)
            x_hat <- self$dec(z)
            list(z = z, x_hat = x_hat)
        }
    )()
    model$to(device = dev)
    opt <- torch::optim_adam(model$parameters, lr = lr)

    n <- nrow(X)
    idx_all <- seq_len(n)
    for (ep in seq_len(epochs)) {
        idx <- sample(idx_all, n)
        batches <- split(idx, ceiling(seq_along(idx) / batch_size))
        for (b in batches) {
            opt$zero_grad()
            out <- model(X_tensor[b, ])
            loss <- torch::nnf_mse_loss(out$x_hat, X_tensor[b, ])
            loss$backward()
            opt$step()
        }
    }

    z <- model(X_tensor)$z$to(device = torch::device("cpu"))$detach()$to(dtype = torch::torch_double())
    emb <- as.matrix(z)
    rownames(emb) <- rownames(X)
    colnames(emb) <- paste0("V", seq_len(dim))
    list(embedding = emb, model = model, method = method)
}

#' Fit diagonal Gaussian embeddings for gene sets from member genes
#'
#' Compute mean and diagonal variance of member gene embeddings for each set.
#'
#' @param gene_embedding Numeric matrix of gene embeddings (genes in rows).
#' @param gene_sets Named list of character vectors (gene IDs).
#' @param eps Lower bound for variances; also used when a set has only one member.
#'
#' @return A list with `mu` and `var` matrices (sets in rows).
#' @export
gsemb_fit_set_gaussians_from_members <- function(gene_embedding,
                                               gene_sets,
                                               eps = 1e-8) {
    gene_sets <- validate_gene_sets(gene_sets)
    E <- as_numeric_matrix(gene_embedding)
    if (is.null(rownames(E))) stop("gene_embedding must have rownames")
    d <- ncol(E)

    set_ids <- names(gene_sets)
    mu <- matrix(NA_real_, length(set_ids), d)
    var <- matrix(NA_real_, length(set_ids), d)
    rownames(mu) <- set_ids
    rownames(var) <- set_ids
    colnames(mu) <- colnames(E)
    colnames(var) <- colnames(E)

    for (sid in set_ids) {
        genes <- intersect(gene_sets[[sid]], rownames(E))
        if (length(genes) == 0) next
        X <- E[genes, , drop = FALSE]
        mu[sid, ] <- colMeans(X)
        if (nrow(X) == 1) {
            var[sid, ] <- rep(eps, d)
        } else {
            v <- apply(X, 2, stats::var)
            var[sid, ] <- pmax(v, eps)
        }
    }
    list(mu = mu, var = var)
}
