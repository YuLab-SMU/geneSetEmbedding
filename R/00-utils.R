as_numeric_matrix <- function(x, id_col = NULL) {
    if (is.matrix(x)) {
        if (!is.numeric(x)) stop("matrix must be numeric")
        storage.mode(x) <- "double"
        return(x)
    }
    if (inherits(x, "Matrix")) {
        mat <- as.matrix(x)
        if (!is.numeric(mat)) stop("matrix must be numeric")
        storage.mode(mat) <- "double"
        return(mat)
    }
    if (is.data.frame(x)) {
        if (is.null(id_col)) {
            if (is.null(rownames(x))) stop("data.frame must have rownames or id_col")
            mat <- as.matrix(x)
        } else {
            if (!id_col %in% names(x)) stop("id_col not found")
            rn <- as.character(x[[id_col]])
            mat <- as.matrix(x[setdiff(names(x), id_col)])
            rownames(mat) <- rn
        }
        storage.mode(mat) <- "double"
        return(mat)
    }
    stop("expected matrix or data.frame")
}

validate_gene_sets <- function(gene_sets) {
    if (!is.list(gene_sets) || is.null(names(gene_sets))) stop("gene_sets must be a named list")
    out <- lapply(gene_sets, function(x) {
        if (!is.character(x)) stop("each gene set must be a character vector")
        unique(x)
    })
    out
}

require_torch <- function() {
    if (!requireNamespace("torch", quietly = TRUE)) stop("torch package is required for this function")
    TRUE
}

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("self"))
}
