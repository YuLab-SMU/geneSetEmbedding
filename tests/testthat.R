library(testthat)
library(geneSetEmbedding)

path <- if (dir.exists("testthat")) {
    "testthat"
} else if (dir.exists(file.path("tests", "testthat"))) {
    file.path("tests", "testthat")
} else {
    system.file("tests/testthat", package = "geneSetEmbedding")
}
test_dir(path)
