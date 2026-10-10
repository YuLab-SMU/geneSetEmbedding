# Soft-membership weighted enrichment (Weighted-ORA / Weighted-GSEA)
# and optional hub-degree bias correction.
#
# Statistic (shared by ORA and GSEA):
#   T_S = sum_g w_gS * stats_g
# where w_gS is a softmax soft membership of gene g to set S
# (optionally reweighted by degree(g)^(-beta) and renormalized).
#
# - Weighted-GSEA: stats_g is continuous (t / logFC / ...).
# - Weighted-ORA:  stats_g is binary 0/1 (DEG indicator).
# - Hub correction: w'_gS ∝ w_gS / degree(g)^beta, then column-normalize.
#   beta = 0 leaves weights unchanged (backward compatible).
#
# gsemb_embedding_enrichment() remains the legacy continuous entry point
# and delegates to the shared driver with restrict_to_members = FALSE
# and degree_beta = 0 by default.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Prefer a precomputed permutation matrix; otherwise try Rcpp; else pure R.
.gsemb_enrichment_null_scores <- function(W, stats_vec, nperm, seed, perm_mat = NULL) {
  stats_vec <- as.numeric(stats_vec)

  if (!is.null(perm_mat)) {
    return(t(crossprod(W, perm_mat)))
  }

  if (.native_routine_available("_geneSetEmbedding_rcpp_enrichment_permutations")) {
    return(rcpp_enrichment_permutations(W, stats_vec, nperm, seed))
  }

  set.seed(seed)
  perm_inner <- vapply(
    seq_len(nperm),
    function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
    FUN.VALUE = numeric(length(stats_vec))
  )

  t(crossprod(W, perm_inner))
}

# Named degree vector aligned to `genes` (missing nodes -> 0).
.gsemb_node_degree_from_adj <- function(adj, genes) {
  if (is.null(adj)) {
    stop("adjacency matrix is required for hub correction (fit$adj or adj=)")
  }
  if (is.null(rownames(adj))) {
    stop("adjacency matrix must have rownames (node / gene IDs)")
  }
  if (inherits(adj, "Matrix")) {
    deg_all <- Matrix::rowSums(adj)
  } else {
    deg_all <- rowSums(as.matrix(adj))
  }
  names(deg_all) <- rownames(adj)
  deg <- as.numeric(deg_all[genes])
  names(deg) <- genes
  deg[is.na(deg)] <- 0
  deg
}

# Softmax gene-to-set scores to soft membership weights W (genes x sets).
# If restrict_to_members=TRUE, non-members get -Inf scores before softmax
# so they do not participate in that set's soft membership.
.gsemb_soft_membership_W <- function(gene_emb,
                                     mu,
                                     var,
                                     gene_sets = NULL,
                                     restrict_to_members = FALSE,
                                     score = c("loglik", "neg_mahalanobis"),
                                     temperature = 1.0,
                                     eps = 1e-8) {
  score <- match.arg(score)
  if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) {
    stop("temperature must be a positive scalar")
  }
  if (isTRUE(restrict_to_members)) {
    if (is.null(gene_sets)) {
      stop("restrict_to_members=TRUE requires gene_sets")
    }
    gene_sets <- validate_gene_sets(gene_sets)
  }

  S <- gsemb_gene_to_set_score(
    gene_embedding = gene_emb,
    set_mu = mu,
    set_var = var,
    score = score,
    eps = eps
  )
  S <- as.matrix(S)
  genes <- rownames(gene_emb)
  set_ids <- rownames(mu)

  if (isTRUE(restrict_to_members)) {
    for (j in seq_along(set_ids)) {
      sid <- set_ids[j]
      if (!sid %in% names(gene_sets)) {
        S[, j] <- -Inf
        next
      }
      members <- intersect(gene_sets[[sid]], genes)
      if (length(members) == 0L) {
        S[, j] <- -Inf
        next
      }
      keep <- genes %in% members
      S[!keep, j] <- -Inf
    }
  }

  col_max <- apply(S, 2, function(col) {
    finite <- col[is.finite(col)]
    if (length(finite) == 0L) {
      return(0)
    }
    max(finite)
  })
  S <- S - rep(col_max, each = nrow(S))
  exp_S <- exp(S / temperature)
  exp_S[!is.finite(exp_S)] <- 0
  col_sum <- colSums(exp_S)
  W <- exp_S / rep(pmax(col_sum, eps), each = nrow(exp_S))
  W[, col_sum <= 0] <- 0
  W[!is.finite(W)] <- 0
  rownames(W) <- genes
  colnames(W) <- set_ids
  W
}

# Hub bias correction: W <- W / degree^beta, then column-normalize.
# beta = 0 returns W unchanged (exact identity, no float drift).
.gsemb_apply_degree_beta <- function(W, degree, beta = 0, eps = 1e-8) {
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta < 0) {
    stop("degree_beta must be a finite non-negative scalar")
  }
  W <- as.matrix(W)
  genes <- rownames(W)
  if (is.null(genes)) {
    stop("W must have rownames")
  }
  if (is.null(names(degree))) {
    if (length(degree) != nrow(W)) {
      stop("degree length must equal nrow(W)")
    }
  } else {
    degree <- as.numeric(degree[genes])
    names(degree) <- genes
  }
  if (length(degree) != nrow(W)) {
    stop("degree and W have incompatible dimensions")
  }

  if (isTRUE(all.equal(beta, 0))) {
    return(W)
  }

  deg_safe <- pmax(as.numeric(degree), 1)
  W2 <- W / (deg_safe^beta)
  W2[!is.finite(W2)] <- 0
  cs <- colSums(W2)
  W2 <- W2 / rep(pmax(cs, eps), each = nrow(W2))
  W2[, cs <= 0] <- 0
  W2[!is.finite(W2)] <- 0
  rownames(W2) <- genes
  colnames(W2) <- colnames(W)
  W2
}

# One batch of sets: soft weights -> optional hub correction -> ES + permutation p.
.gsemb_weighted_enrichment_sets <- function(gene_stats,
                                            gene_emb,
                                            genes,
                                            mu,
                                            var,
                                            gene_sets,
                                            restrict_to_members,
                                            score,
                                            temperature,
                                            nperm,
                                            alternative,
                                            seed,
                                            eps,
                                            top_genes,
                                            perm_mat,
                                            degree = NULL,
                                            degree_beta = 0) {
  W <- .gsemb_soft_membership_W(
    gene_emb = gene_emb,
    mu = mu,
    var = var,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    eps = eps
  )
  if (!is.null(degree) || !isTRUE(all.equal(degree_beta, 0))) {
    if (is.null(degree) && !isTRUE(all.equal(degree_beta, 0))) {
      stop("degree_beta > 0 requires a degree vector")
    }
    if (!is.null(degree)) {
      W <- .gsemb_apply_degree_beta(W, degree = degree, beta = degree_beta, eps = eps)
    }
  }

  # gene_stats must already be ordered as rownames(gene_emb)
  stats_vec <- as.numeric(gene_stats)
  es <- as.numeric(crossprod(W, stats_vec))
  names(es) <- colnames(W)
  n_es <- length(es)
  set_ids <- colnames(W)

  set_size <- rep(NA_integer_, n_es)
  if (!is.null(gene_sets)) {
    set_size <- vapply(set_ids, function(sid) {
      if (!sid %in% names(gene_sets)) {
        return(NA_integer_)
      }
      length(intersect(gene_sets[[sid]], genes))
    }, integer(1))
  }

  # Distinguish untestable sets from non-significant ones (gap doc §2.3).
  # Priority: params > coverage > weights > ok.
  w_col_sum <- colSums(W)
  status <- character(n_es)
  for (j in seq_len(n_es)) {
    if (!all(is.finite(mu[j, ])) || !all(is.finite(var[j, ]))) {
      status[j] <- "untestable_params"
    } else if (!is.null(gene_sets) && !is.na(set_size[j]) && set_size[j] == 0L) {
      status[j] <- "untestable_coverage"
    } else if (!is.finite(w_col_sum[j]) || w_col_sum[j] <= 0) {
      status[j] <- "untestable_weights"
    } else {
      status[j] <- "ok"
    }
  }
  testable <- status == "ok"

  z <- rep(NA_real_, n_es)
  pvals <- rep(NA_real_, n_es)

  if (nperm > 0L) {
    null_scores <- .gsemb_enrichment_null_scores(W, stats_vec, nperm, seed, perm_mat)
    null_mean <- colMeans(null_scores)
    null_sd <- pmax(apply(null_scores, 2, stats::sd), eps)
    z <- (es - null_mean) / null_sd

    es_mat <- matrix(es, nrow = nperm, ncol = n_es, byrow = TRUE)
    if (alternative == "greater") {
      k <- colSums(null_scores >= es_mat)
    } else if (alternative == "less") {
      k <- colSums(null_scores <= es_mat)
    } else {
      cen_mat <- matrix(null_mean, nrow = nperm, ncol = n_es, byrow = TRUE)
      shift_mat <- matrix(abs(es - null_mean), nrow = nperm, ncol = n_es, byrow = TRUE)
      k <- colSums(abs(null_scores - cen_mat) >= shift_mat)
    }
    # Application / Phipson-Smyth style Monte Carlo p-value
    pvals <- (k + 1) / (nperm + 1)
  }

  # Untestable sets must not look like ES=0 / non-significant.
  es[!testable] <- NA_real_
  z[!testable] <- NA_real_
  pvals[!testable] <- NA_real_

  core <- rep(NA_character_, n_es)
  if (top_genes > 0L) {
    for (j in seq_len(n_es)) {
      if (!testable[j]) {
        next
      }
      wj <- W[, j]
      pos <- which(wj > 0)
      if (length(pos) == 0L) {
        next
      }
      ord <- pos[order(wj[pos], decreasing = TRUE)]
      ord <- ord[seq_len(min(top_genes, length(ord)))]
      core[j] <- paste0(rownames(W)[ord], collapse = "/")
    }
  }

  data.frame(
    ID = set_ids,
    ES = as.numeric(es),
    z = as.numeric(z),
    pvalue = as.numeric(pvals),
    p.adjust = NA_real_,
    setSize = as.integer(set_size),
    core_enrichment = core,
    degree_beta = as.numeric(degree_beta),
    status = status,
    stringsAsFactors = FALSE
  )
}

# Shared driver for Weighted-ORA / Weighted-GSEA / legacy enrichment.
.gsemb_weighted_enrichment_driver <- function(gene_stats,
                                              x,
                                              sets = NULL,
                                              gene_sets = NULL,
                                              restrict_to_members = FALSE,
                                              score = c("loglik", "neg_mahalanobis"),
                                              temperature = 1.0,
                                              nperm = 1000,
                                              alternative = c("two.sided", "greater", "less"),
                                              seed = 1,
                                              eps = 1e-8,
                                              top_genes = 30,
                                              degree_beta = 0,
                                              adj = NULL) {
  score <- match.arg(score)
  alternative <- match.arg(alternative)
  if (!inherits(x, "gsemb_embedding")) {
    stop("x must be a gsemb_embedding object")
  }
  if (!is.numeric(gene_stats) || is.null(names(gene_stats))) {
    stop("gene_stats must be a named numeric vector")
  }
  if (!is.numeric(temperature) || length(temperature) != 1 || temperature <= 0) {
    stop("temperature must be a positive scalar")
  }
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) {
    stop("nperm must be a non-negative integer")
  }
  nperm <- as.integer(nperm)
  if (!is.numeric(top_genes) || length(top_genes) != 1 || top_genes < 0) {
    stop("top_genes must be a non-negative integer")
  }
  top_genes <- as.integer(top_genes)
  if (!is.numeric(degree_beta) || length(degree_beta) != 1L ||
        !is.finite(degree_beta) || degree_beta < 0) {
    stop("degree_beta must be a finite non-negative scalar")
  }
  if (isTRUE(restrict_to_members) && is.null(gene_sets)) {
    stop("restrict_to_members=TRUE requires gene_sets")
  }

  gene_emb <- x$gene_embedding
  if (is.null(rownames(gene_emb))) {
    stop("gene embedding must have rownames")
  }
  mu <- as_numeric_matrix(x$set_mu)
  var <- as_numeric_matrix(x$set_var)
  if (is.null(rownames(mu)) || is.null(rownames(var))) {
    stop("set embedding must have rownames")
  }
  if (nrow(mu) != nrow(var)) {
    stop("set_mu and set_var must have the same number of rows")
  }
  if (length(rownames(mu)) != nrow(mu)) {
    stop(
      "set_mu rownames length (", length(rownames(mu)),
      ") != nrow (", nrow(mu), ")"
    )
  }

  genes <- intersect(names(gene_stats), rownames(gene_emb))
  if (length(genes) < 2) {
    stop("not enough genes overlap between gene_stats and gene_embedding")
  }
  gene_emb <- gene_emb[genes, , drop = FALSE]
  # Align stats to embedding row order (same as legacy enrichment).
  gene_stats <- gene_stats[rownames(gene_emb)]
  storage.mode(gene_stats) <- "double"
  genes <- rownames(gene_emb)

  if (is.null(adj)) {
    adj <- x$adj
  }
  degree <- NULL
  if (!is.null(adj) || !isTRUE(all.equal(degree_beta, 0))) {
    degree <- .gsemb_node_degree_from_adj(adj, genes)
  }

  set_names <- rownames(mu)
  if (is.null(sets)) {
    set_idx <- seq_len(nrow(mu))
  } else {
    set_idx <- match(as.character(sets), set_names, nomatch = 0L)
    if (any(set_idx == 0L)) {
      bad <- sets[set_idx == 0L]
      stop(
        "some sets are not in set_mu rownames; first missing: ",
        paste(utils::head(bad, 3L), collapse = ", ")
      )
    }
  }
  if (length(set_idx) == 0L) {
    stop("no sets found in embedding")
  }

  # Subset first, then batch with local 1..n indices (avoids indexing bugs
  # when `sets` is a non-contiguous subset of rownames).
  mu <- mu[set_idx, , drop = FALSE]
  var <- var[set_idx, , drop = FALSE]
  n_local <- nrow(mu)
  local_idx <- seq_len(n_local)

  if (!is.null(gene_sets)) {
    gene_sets <- validate_gene_sets(gene_sets)
  }

  # One shared perm_mat for all batches so p-values stay coherent.
  perm_mat <- NULL
  if (nperm > 0L) {
    stats_vec <- as.numeric(gene_stats)
    set.seed(seed)
    perm_mat <- vapply(
      seq_len(nperm),
      function(b) sample(stats_vec, length(stats_vec), replace = FALSE),
      FUN.VALUE = numeric(length(stats_vec))
    )
  }

  # Auto-batch when n_genes * n_sets is large (same threshold as legacy code).
  n_genes <- nrow(gene_emb)
  sets_per_batch <- if (as.numeric(n_genes) * as.numeric(n_local) > 2e7) {
    max(50L, as.integer(floor(2e7 / max(n_genes, 1L))))
  } else {
    n_local
  }
  if (!is.finite(sets_per_batch) || sets_per_batch < 1L) {
    sets_per_batch <- n_local
  }

  idx_batches <- split(local_idx, ceiling(seq_along(local_idx) / sets_per_batch))
  parts <- lapply(idx_batches, function(batch_idx) {
    .gsemb_weighted_enrichment_sets(
      gene_stats = gene_stats,
      gene_emb = gene_emb,
      genes = genes,
      mu = mu[batch_idx, , drop = FALSE],
      var = var[batch_idx, , drop = FALSE],
      gene_sets = gene_sets,
      restrict_to_members = restrict_to_members,
      score = score,
      temperature = temperature,
      nperm = nperm,
      alternative = alternative,
      seed = seed,
      eps = eps,
      top_genes = top_genes,
      perm_mat = perm_mat,
      degree = degree,
      degree_beta = degree_beta
    )
  })

  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out$p.adjust <- NA_real_
  ok <- !is.null(out$status) & out$status == "ok" & is.finite(out$pvalue)
  if (any(ok)) {
    out$p.adjust[ok] <- stats::p.adjust(out$pvalue[ok], method = "BH")
  }
  out
}

# ---------------------------------------------------------------------------
# Exported API
# ---------------------------------------------------------------------------

#' Soft-membership enrichment (legacy continuous entry point)
#'
#' Computes set enrichment scores from gene-level statistics using soft
#' membership weights derived from diagonal Gaussian set embeddings.
#' By default this matches the historical behaviour:
#' `restrict_to_members = FALSE` and `degree_beta = 0` (no hub correction).
#' Prefer [gsemb_weighted_gsea()] / [gsemb_weighted_ora()] for the explicit
#' Weighted-GSEA / Weighted-ORA interfaces.
#'
#' @param gene_stats Named numeric vector of gene-level statistics.
#' @param x A `gsemb_embedding` object.
#' @param sets Optional character vector of set IDs to test.
#' @param gene_sets Named list of gene-set members (set-size reporting; required
#'   when `restrict_to_members = TRUE`).
#' @param score `"loglik"` or `"neg_mahalanobis"`.
#' @param temperature Softmax temperature for soft membership weights.
#' @param nperm Number of permutations for the null distribution.
#' @param alternative `"two.sided"`, `"greater"`, or `"less"`.
#' @param seed Random seed for permutations.
#' @param eps Numerical floor for variances / empty columns.
#' @param top_genes Number of top positive-weight genes reported in
#'   `core_enrichment` (not a classical GSEA leading edge).
#' @param restrict_to_members If `TRUE`, only original set members enter that
#'   set's soft membership (non-members get weight 0). Default `FALSE`.
#' @param degree_beta Hub penalty exponent for `w / degree^beta`
#'   (default `0` = no correction). Typical probes: `0, 0.25, 0.5, 1`.
#' @param adj Optional adjacency for degrees; defaults to `x$adj`.
#'
#' @return A data.frame with `ID`, `ES`, `z`, `pvalue`, `p.adjust`, `setSize`,
#'   `core_enrichment`, `degree_beta`, and `status`. Untestable sets
#'   (`untestable_params`, `untestable_coverage`, or `untestable_weights`) have
#'   `ES`/`z`/`pvalue`/`p.adjust`/`core_enrichment` as `NA` and are excluded
#'   from BH; only `status == "ok"` rows enter `p.adjust`.
#' @seealso [gsemb_weighted_gsea()], [gsemb_weighted_ora()], [gsemb_sweep_hub_beta()]
#' @examples
#' \dontrun{
#' enrich <- gsemb_embedding_enrichment(
#'   gene_stats = stats,
#'   x = fit,
#'   gene_sets = gene_sets,
#'   nperm = 100,
#'   alternative = "greater"
#' )
#' }
#' @export
gsemb_embedding_enrichment <- function(gene_stats,
                                       x,
                                       sets = NULL,
                                       gene_sets = NULL,
                                       score = c("loglik", "neg_mahalanobis"),
                                       temperature = 1.0,
                                       nperm = 1000,
                                       alternative = c("two.sided", "greater", "less"),
                                       seed = 1,
                                       eps = 1e-8,
                                       top_genes = 30,
                                       restrict_to_members = FALSE,
                                       degree_beta = 0,
                                       adj = NULL) {
  .gsemb_weighted_enrichment_driver(
    gene_stats = gene_stats,
    x = x,
    sets = sets,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    nperm = nperm,
    alternative = alternative,
    seed = seed,
    eps = eps,
    top_genes = top_genes,
    degree_beta = degree_beta,
    adj = adj
  )
}

#' Weighted-GSEA with soft membership weights
#'
#' Continuous gene statistics \code{r_g} (e.g. t-statistic or logFC) yield
#' \code{T_S = sum_g w_gS * r_g}. Optional hub correction reweights by
#' \code{w_gS / degree(g)^beta}.
#'
#' @inheritParams gsemb_embedding_enrichment
#' @param r Named numeric vector of continuous gene-level statistics.
#' @return Same columns as [gsemb_embedding_enrichment()], including `status`;
#'   BH is applied only to testable (`status == "ok"`) rows.
#' @export
gsemb_weighted_gsea <- function(r,
                                x,
                                sets = NULL,
                                gene_sets = NULL,
                                score = c("loglik", "neg_mahalanobis"),
                                temperature = 1.0,
                                nperm = 1000,
                                alternative = c("two.sided", "greater", "less"),
                                seed = 1,
                                eps = 1e-8,
                                top_genes = 30,
                                restrict_to_members = FALSE,
                                degree_beta = 0,
                                adj = NULL) {
  if (!is.numeric(r) || is.null(names(r))) {
    stop("r must be a named numeric vector of continuous gene statistics")
  }
  .gsemb_weighted_enrichment_driver(
    gene_stats = r,
    x = x,
    sets = sets,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    nperm = nperm,
    alternative = alternative,
    seed = seed,
    eps = eps,
    top_genes = top_genes,
    degree_beta = degree_beta,
    adj = adj
  )
}

#' Weighted-ORA with soft membership weights
#'
#' Binary labels \code{y_g} in \code{c(0, 1)} (e.g. differential-expression
#' indicators) yield \code{T_S = sum_g w_gS * y_g}. Same soft-weight / hub /
#' permutation machinery as [gsemb_weighted_gsea()].
#'
#' @inheritParams gsemb_embedding_enrichment
#' @param y Named numeric vector with values in \code{c(0, 1)} only.
#' @return Same columns as [gsemb_embedding_enrichment()], including `status`;
#'   BH is applied only to testable (`status == "ok"`) rows.
#' @export
gsemb_weighted_ora <- function(y,
                               x,
                               sets = NULL,
                               gene_sets = NULL,
                               score = c("loglik", "neg_mahalanobis"),
                               temperature = 1.0,
                               nperm = 1000,
                               alternative = c("two.sided", "greater", "less"),
                               seed = 1,
                               eps = 1e-8,
                               top_genes = 30,
                               restrict_to_members = FALSE,
                               degree_beta = 0,
                               adj = NULL) {
  if (!is.numeric(y) || is.null(names(y))) {
    stop("y must be a named numeric vector of 0/1 labels")
  }
  y_names <- names(y)
  y_num <- as.numeric(y)
  names(y_num) <- y_names
  uniq <- unique(y_num[is.finite(y_num)])
  if (length(uniq) == 0L || !all(uniq %in% c(0, 1))) {
    stop("Weighted-ORA requires y to contain only 0 and 1")
  }
  storage.mode(y_num) <- "double"

  .gsemb_weighted_enrichment_driver(
    gene_stats = y_num,
    x = x,
    sets = sets,
    gene_sets = gene_sets,
    restrict_to_members = restrict_to_members,
    score = score,
    temperature = temperature,
    nperm = nperm,
    alternative = alternative,
    seed = seed,
    eps = eps,
    top_genes = top_genes,
    degree_beta = degree_beta,
    adj = adj
  )
}

#' Sweep hub-correction strengths (degree_beta)
#'
#' Runs Weighted-GSEA or Weighted-ORA for each beta in `betas` (default
#' `c(0, 0.25, 0.5, 1)` as commonly used in grant-style probes; any
#' non-negative values are allowed).
#'
#' @param gene_stats Named numeric vector; must be 0/1 when `mode = "ora"`.
#' @param x A `gsemb_embedding` object.
#' @param betas Numeric vector of non-negative hub penalties.
#' @param mode `"gsea"` or `"ora"`.
#' @param ... Passed to [gsemb_weighted_gsea()] or [gsemb_weighted_ora()]
#'   (e.g. `gene_sets`, `nperm`, `seed`, `restrict_to_members`).
#'
#' @return A list with:
#'   \describe{
#'     \item{results}{Named list of per-beta data.frames}
#'     \item{combined}{Row-bound long table}
#'     \item{compare}{Wide ES table plus `deltaES_vs0_*` columns}
#'   }
#' @export
gsemb_sweep_hub_beta <- function(gene_stats,
                                 x,
                                 betas = c(0, 0.25, 0.5, 1),
                                 mode = c("gsea", "ora"),
                                 ...) {
  mode <- match.arg(mode)
  betas <- as.numeric(betas)
  if (any(!is.finite(betas) | betas < 0)) {
    stop("betas must be finite and non-negative")
  }

  fun <- if (mode == "ora") gsemb_weighted_ora else gsemb_weighted_gsea
  results <- lapply(betas, function(b) {
    fun(gene_stats, x, degree_beta = b, ...)
  })
  names(results) <- paste0("beta_", betas)

  combined <- do.call(rbind, results)
  rownames(combined) <- NULL

  es_list <- lapply(seq_along(betas), function(i) {
    d <- results[[i]]
    v <- d$ES
    names(v) <- d$ID
    v
  })
  all_ids <- unique(unlist(lapply(es_list, names), use.names = FALSE))
  compare <- data.frame(ID = all_ids, stringsAsFactors = FALSE)
  for (i in seq_along(betas)) {
    compare[[paste0("ES_beta_", betas[i])]] <- as.numeric(es_list[[i]][all_ids])
  }
  if ("ES_beta_0" %in% names(compare)) {
    for (b in betas[betas != 0]) {
      col <- paste0("ES_beta_", b)
      compare[[paste0("deltaES_vs0_", b)]] <- compare[[col]] - compare$ES_beta_0
    }
  }

  list(results = results, combined = combined, compare = compare)
}
