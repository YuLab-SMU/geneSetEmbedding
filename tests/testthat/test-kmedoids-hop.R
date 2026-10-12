library(testthat)
library(geneSetEmbedding)

test_that("kmedoids medoids match the pre-reuse landmarks", {
  path <- data.frame(
    node1 = c("A", "B", "C", "D", "E", "F"),
    node2 = c("B", "C", "D", "E", "F", "G")
  )
  adj <- gsemb_build_graph(path, directed = FALSE)
  lm <- gsemb_select_landmarks(
    adj,
    k = 2,
    method = "kmedoids",
    seed = 1,
    kmedoids_m = 50,
    kmedoids_max_iter = 10
  )
  expect_identical(lm, c("C", "F"))
})

test_that("kmedoids medoids match the pre-reuse weighted landmarks", {
  edges <- data.frame(
    node1 = c("A", "B", "C", "D", "A", "B"),
    node2 = c("B", "C", "D", "E", "E", "D"),
    weight = c(100, 100, 100, 100, 1, 1)
  )
  adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)
  lm <- gsemb_select_landmarks(
    adj,
    k = 2,
    method = "kmedoids",
    seed = 1,
    kmedoids_m = 50,
    kmedoids_max_iter = 10
  )
  expect_identical(lm, c("D", "A"))
})

test_that("k at least n returns every node of the largest component", {
  edges <- data.frame(
    node1 = c("A", "B", "C", "X"),
    node2 = c("B", "C", "D", "Y")
  )
  adj <- gsemb_build_graph(edges, directed = FALSE)
  lm <- gsemb_select_landmarks(
    adj,
    k = 10,
    method = "kmedoids",
    seed = 1,
    kmedoids_m = 50,
    kmedoids_max_iter = 10
  )
  expect_identical(lm, c("A", "B", "C", "D"))
})
