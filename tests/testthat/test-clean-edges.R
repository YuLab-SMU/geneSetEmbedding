library(testthat)
library(geneSetEmbedding)

test_that("gsemb_clean_edges drops self-loops and collapses undirected pairs", {
  edges <- data.frame(
    node1 = c("B", "A", "A", "C"),
    node2 = c("A", "B", "A", "D"),
    weight = c(800, 900, 50, 500),
    stringsAsFactors = FALSE
  )
  out <- gsemb_clean_edges(edges, weight = "weight", score_cutoff = 700)

  expect_equal(nrow(out), 1L)
  expect_identical(out$node1, "A")
  expect_identical(out$node2, "B")
  expect_equal(out$weight, 900)
})

test_that("gsemb_clean_edges keeps a one-way edge after endpoint ordering", {
  edges <- data.frame(
    node1 = "B",
    node2 = "A",
    weight = 10,
    stringsAsFactors = FALSE
  )
  out <- gsemb_clean_edges(edges, weight = "weight")

  expect_equal(nrow(out), 1L)
  expect_identical(out$node1, "A")
  expect_identical(out$node2, "B")
})

test_that("gsemb_clean_edges works without a weight column", {
  edges <- data.frame(
    node1 = c("A", "B"),
    node2 = c("B", "A"),
    stringsAsFactors = FALSE
  )
  out <- gsemb_clean_edges(edges)

  expect_equal(nrow(out), 1L)
  expect_identical(out$node1, "A")
  expect_identical(out$node2, "B")
})

test_that("score_cutoff requires weight", {
  edges <- data.frame(node1 = "A", node2 = "B")
  expect_error(
    gsemb_clean_edges(edges, score_cutoff = 700),
    "weight must be set"
  )
})

test_that("missing edge endpoints are counted and skipped", {
  edges <- data.frame(
    node1 = c("A", NA, "B", ""),
    node2 = c("B", "C", NA, "D"),
    weight = c(800, 900, 850, 750)
  )
  expect_warning(
    out <- gsemb_clean_edges(edges, weight = "weight"),
    "Dropped 3 edges with missing endpoints"
  )
  expect_equal(nrow(out), 1L)
  expect_identical(out$node1, "A")
  expect_identical(out$node2, "B")
})
