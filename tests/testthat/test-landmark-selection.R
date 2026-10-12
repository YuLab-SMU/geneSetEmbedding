library(testthat)
library(geneSetEmbedding)

test_that("gsemb_fit forwards landmark_method into landmark selection", {
  edges <- data.frame(
    node1 = c("A", "B", "C", "D"),
    node2 = c("B", "C", "D", "E"),
    weight = c(1, 1, 1, 1)
  )
  gene_sets <- list(S = c("A", "B", "C"))
  fit_one <- function(landmark_method) {
    gsemb_fit(
      edges,
      gene_sets,
      weight = "weight",
      method = "svd",
      dim = 1,
      k = 1,
      max_iter = 10,
      landmark_method = landmark_method,
      betweenness_cutoff = -1,
      kmedoids_m = 10
    )
  }

  expect_identical(fit_one("degree")$landmarks, "B")
  expect_identical(fit_one("degree")$landmark_method, "degree")
  expect_identical(fit_one("betweenness")$landmarks, "C")
  expect_identical(fit_one("kmedoids")$landmarks, "C")
  expect_error(fit_one("weighted_degree"), "should be one of")
})

test_that("degree selects nodes by incident edge-weight sum", {
  edges <- data.frame(
    node1 = c("A", "A", "D", "D", "D"),
    node2 = c("B", "C", "B", "C", "E"),
    weight = c(10, 10, 3, 3, 3)
  )
  adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)

  expect_identical(
    gsemb_select_landmarks(adj, k = 1, method = "degree"),
    "A"
  )
})

test_that("path-based landmark methods reject nonpositive and nonfinite weights", {
  edges <- data.frame(
    node1 = c("A", "B"), node2 = c("B", "C"), weight = c(800, -1)
  )
  adj <- gsemb_build_graph(edges, weight = "weight")
  for (method in c("betweenness", "kmedoids")) {
    expect_error(
      gsemb_select_landmarks(adj, k = 1, method = method),
      "finite and positive"
    )
  }

  adj["B", "C"] <- Inf
  adj["C", "B"] <- Inf
  expect_error(
    gsemb_select_landmarks(adj, k = 1, method = "betweenness"),
    "finite and positive"
  )

  positive <- gsemb_build_graph(
    transform(edges, weight = c(800, 2000)), weight = "weight"
  )
  for (method in c("betweenness", "kmedoids")) {
    expect_length(gsemb_select_landmarks(positive, k = 1, method = method), 1L)
  }
})

test_that("betweenness prefers bridge nodes on a path graph", {
  skip_if_not_installed("igraph")
  # A-B-C-D-E: middle node C has highest betweenness
  edges <- data.frame(
    node1 = c("A", "B", "C", "D"),
    node2 = c("B", "C", "D", "E")
  )
  adj <- gsemb_build_graph(edges, directed = FALSE)

  expect_identical(
    gsemb_select_landmarks(adj, k = 1, method = "betweenness"),
    "C"
  )
  expect_identical(
    gsemb_select_landmarks(adj, k = 1, method = "betweenness", betweenness_cutoff = -1),
    "C"
  )
})

test_that("betweenness prefers the high-weight path bridge", {
  skip_if_not_installed("igraph")
  # Two A--*--C paths: high weights via B, low via D → weighted picks B
  edges <- data.frame(
    node1 = c("A", "B", "A", "D"),
    node2 = c("B", "C", "D", "C"),
    weight = c(100, 100, 1, 1)
  )
  adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)

  expect_identical(
    gsemb_select_landmarks(
      adj,
      k = 1,
      method = "betweenness",
      betweenness_cutoff = -1
    ),
    "B"
  )
})

test_that("kmedoids on a path with k=1 converges to the center", {
  skip_if_not_installed("igraph")
  edges <- data.frame(
    node1 = c("A", "B", "C", "D"),
    node2 = c("B", "C", "D", "E")
  )
  adj <- gsemb_build_graph(edges, directed = FALSE)

  expect_identical(
    gsemb_select_landmarks(adj, k = 1, method = "kmedoids", kmedoids_m = 10),
    "C"
  )
})

test_that("kmedoids k=2 covers two cliques in one component", {
  skip_if_not_installed("igraph")
  # Two K4s joined by a path (one CC). Disconnected equal cliques would
  # collapse to the first largest CC under the current helper.
  clique_edges <- function(ids) {
    data.frame(
      node1 = c(ids[1], ids[1], ids[1], ids[2], ids[2], ids[3]),
      node2 = c(ids[2], ids[3], ids[4], ids[3], ids[4], ids[4])
    )
  }
  a <- c("A1", "A2", "A3", "A4")
  b <- c("B1", "B2", "B3", "B4")
  edges <- rbind(
    clique_edges(a),
    data.frame(node1 = c("A1", "X", "Y"), node2 = c("X", "Y", "B1")),
    clique_edges(b)
  )
  adj <- gsemb_build_graph(edges, directed = FALSE)
  lm <- gsemb_select_landmarks(
    adj,
    k = 2,
    method = "kmedoids",
    seed = 1,
    kmedoids_m = 20
  )

  expect_true(any(lm %in% a) && any(lm %in% b))
})

test_that("kmedoids is reproducible and returns k unique node IDs", {
  skip_if_not_installed("igraph")
  edges <- data.frame(
    node1 = c("A", "B", "C", "D"),
    node2 = c("B", "C", "D", "E")
  )
  adj <- gsemb_build_graph(edges, directed = FALSE)
  a <- gsemb_select_landmarks(adj, k = 2, method = "kmedoids", seed = 7)
  b <- gsemb_select_landmarks(adj, k = 2, method = "kmedoids", seed = 7)

  expect_identical(a, b)
  expect_length(a, 2L)
  expect_equal(length(unique(a)), 2L)
  expect_true(all(a %in% rownames(adj)))
})

test_that("kmedoids larger than the largest CC returns that component", {
  skip_if_not_installed("igraph")
  edges <- data.frame(
    node1 = c("A", "B", "C"),
    node2 = c("B", "C", "A")
  )
  adj <- gsemb_build_graph(
    edges,
    nodes = c("A", "B", "C", "X"),
    directed = FALSE
  )
  lm <- gsemb_select_landmarks(adj, k = 10, method = "kmedoids")

  expect_setequal(lm, c("A", "B", "C"))
})

test_that("kmedoids prefers endpoints of the high-weight edge", {
  skip_if_not_installed("igraph")
  # Triangle: A-C is heavy. Hop costs are tied; weighted optima are A or C (not B).
  edges <- data.frame(
    node1 = c("A", "B", "A"),
    node2 = c("B", "C", "C"),
    weight = c(1, 1, 100)
  )
  adj <- gsemb_build_graph(edges, weight = "weight", directed = FALSE)

  picks <- unique(vapply(1:30, function(s) {
    gsemb_select_landmarks(
      adj, k = 1, method = "kmedoids", seed = s, kmedoids_m = 10
    )
  }, character(1)))
  expect_true(all(picks %in% c("A", "C")))
  expect_false("B" %in% picks)
})
