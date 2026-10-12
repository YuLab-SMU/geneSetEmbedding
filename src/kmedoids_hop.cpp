// K-Medoids landmarks: farthest-first + restricted PAM.
// Edge cost is 1 / weight (Dijkstra). Do NOT add // [[Rcpp::export]];
// register in geneSetEmbedding_init.cpp.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <random>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();

struct Edge {
  int to;
  double cost;
};

void dijkstra_from(int src,
                   const std::vector<std::vector<Edge>>& adj,
                   std::vector<double>& dist) {
  std::fill(dist.begin(), dist.end(), kInf);
  using Node = std::pair<double, int>;
  std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;
  dist[static_cast<std::size_t>(src)] = 0.0;
  pq.push(Node(0.0, src));
  while (!pq.empty()) {
    const double du = pq.top().first;
    const int u = pq.top().second;
    pq.pop();
    if (du > dist[static_cast<std::size_t>(u)]) {
      continue;
    }
    for (const Edge& e : adj[static_cast<std::size_t>(u)]) {
      const double nd = du + e.cost;
      if (nd < dist[static_cast<std::size_t>(e.to)]) {
        dist[static_cast<std::size_t>(e.to)] = nd;
        pq.push(Node(nd, e.to));
      }
    }
  }
}


double total_cost(const std::vector<std::vector<double>>& dist_cols, int n, int k) {
  double cost = 0.0;
  for (int v = 0; v < n; ++v) {
    double best = kInf;
    for (int j = 0; j < k; ++j) {
      best = std::min(best, dist_cols[static_cast<std::size_t>(j)][static_cast<std::size_t>(v)]);
    }
    cost += best;
  }
  return cost;
}

double trial_cost(const std::vector<double>& best_other,
                  const std::vector<double>& dnew,
                  int n) {
  double cost = 0.0;
  for (int v = 0; v < n; ++v) {
    cost += std::min(dnew[static_cast<std::size_t>(v)],
                     best_other[static_cast<std::size_t>(v)]);
  }
  return cost;
}

void assign_clusters(const std::vector<std::vector<double>>& dist_cols,
                     int n,
                     int k,
                     std::vector<int>& assign) {
  assign.resize(static_cast<std::size_t>(n));
  for (int v = 0; v < n; ++v) {
    int best_j = 0;
    double best_d = dist_cols[0][static_cast<std::size_t>(v)];
    for (int j = 1; j < k; ++j) {
      const double d = dist_cols[static_cast<std::size_t>(j)][static_cast<std::size_t>(v)];
      if (d < best_d) {
        best_d = d;
        best_j = j;
      }
    }
    assign[static_cast<std::size_t>(v)] = best_j;
  }
}

}  // namespace

// i, p, x: 0-based CSC of undirected adjacency (column j lists neighbors / weights).
// Edge cost = 1/weight. Returns 0-based medoid vertex indices (length k).
RcppExport SEXP _geneSetEmbedding_kmedoids_hop(SEXP i_,
                                               SEXP p_,
                                               SEXP x_,
                                               SEXP k_,
                                               SEXP seed_,
                                               SEXP m_,
                                               SEXP max_iter_) {
  BEGIN_RCPP
  const IntegerVector i(i_);
  const IntegerVector p(p_);
  const NumericVector x(x_);
  const int k_req = as<int>(k_);
  const int seed = as<int>(seed_);
  const int m = as<int>(m_);
  const int max_iter = as<int>(max_iter_);

  const int n = p.size() - 1;
  if (n <= 0) {
    stop("kmedoids_hop: empty graph");
  }
  if (x.size() != i.size()) {
    stop("kmedoids_hop: x and i length mismatch");
  }
  int k = std::min(k_req, n);
  if (k <= 0) {
    stop("kmedoids_hop: k must be positive");
  }

  std::vector<std::vector<Edge>> adj(static_cast<std::size_t>(n));
  std::vector<double> strength(static_cast<std::size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    const int start = p[j];
    const int end = p[j + 1];
    auto& nbrs = adj[static_cast<std::size_t>(j)];
    nbrs.reserve(static_cast<std::size_t>(end - start));
    for (int e = start; e < end; ++e) {
      const double w = x[e];
      const double cost = 1.0 / std::max(w, 1e-15);
      nbrs.push_back(Edge{i[e], cost});
      strength[static_cast<std::size_t>(j)] += w;
    }
  }

  if (k >= n) {
    IntegerVector out(n);
    for (int v = 0; v < n; ++v) {
      out[v] = v;
    }
    return out;
  }

  std::mt19937 rng(static_cast<std::uint32_t>(seed));
  std::uniform_int_distribution<int> uni(0, n - 1);

  std::vector<std::vector<double>> dist_cols(static_cast<std::size_t>(k),
                                             std::vector<double>(static_cast<std::size_t>(n)));
  std::vector<int> medoids(static_cast<std::size_t>(k), -1);
  std::vector<char> selected(static_cast<std::size_t>(n), 0);
  std::vector<double> min_dist(static_cast<std::size_t>(n), kInf);
  // ponytail: one n-vector per distinct source. First PAM pass still computes new candidates.
  // A full-graph cache is n columns (~2GB at n=16185); only queried sources are stored.
  std::vector<char> have_dist(static_cast<std::size_t>(n), 0);
  std::vector<std::vector<double>> dist_cache(static_cast<std::size_t>(n));
  auto cached_sssp = [&](int src) -> const std::vector<double>& {
    const std::size_t s = static_cast<std::size_t>(src);
    if (!have_dist[s]) {
      dist_cache[s].resize(static_cast<std::size_t>(n));
      dijkstra_from(src, adj, dist_cache[s]);
      have_dist[s] = 1;
    }
    return dist_cache[s];
  };

  medoids[0] = uni(rng);
  selected[static_cast<std::size_t>(medoids[0])] = 1;
  dist_cols[0] = cached_sssp(medoids[0]);
  min_dist = dist_cols[0];

  for (int t = 1; t < k; ++t) {
    int best_v = -1;
    double best_d = -1.0;
    for (int v = 0; v < n; ++v) {
      if (selected[static_cast<std::size_t>(v)]) {
        continue;
      }
      const double d = min_dist[static_cast<std::size_t>(v)];
      if (d > best_d) {
        best_d = d;
        best_v = v;
      }
    }
    if (best_v < 0) {
      for (int v = 0; v < n; ++v) {
        if (!selected[static_cast<std::size_t>(v)]) {
          best_v = v;
          break;
        }
      }
    }
    medoids[static_cast<std::size_t>(t)] = best_v;
    selected[static_cast<std::size_t>(best_v)] = 1;
    dist_cols[static_cast<std::size_t>(t)] = cached_sssp(best_v);
    for (int v = 0; v < n; ++v) {
      min_dist[static_cast<std::size_t>(v)] = std::min(
          min_dist[static_cast<std::size_t>(v)],
          dist_cols[static_cast<std::size_t>(t)][static_cast<std::size_t>(v)]);
    }
  }

  std::vector<int> assign;
  // ponytail: Top-M by strength; raise m for fuller PAM.
  const int m_cap = std::max(1, m);

  for (int iter = 0; iter < max_iter; ++iter) {
    assign_clusters(dist_cols, n, k, assign);
    double cost = total_cost(dist_cols, n, k);
    bool improved = false;

    for (int j = 0; j < k; ++j) {
      std::vector<int> members;
      members.reserve(static_cast<std::size_t>(n / std::max(k, 1)));
      for (int v = 0; v < n; ++v) {
        if (assign[static_cast<std::size_t>(v)] == j) {
          members.push_back(v);
        }
      }
      if (members.size() <= 1) {
        continue;
      }

      std::sort(members.begin(), members.end(), [&](int a, int b) {
        const double sa = strength[static_cast<std::size_t>(a)];
        const double sb = strength[static_cast<std::size_t>(b)];
        if (sa != sb) {
          return sa > sb;
        }
        return a < b;
      });
      const int n_cand = std::min(m_cap, static_cast<int>(members.size()));

      // ponytail: one best_other column per replaced center. Candidate count is still m.
      std::vector<double> best_other(static_cast<std::size_t>(n), kInf);
      for (int t = 0; t < k; ++t) {
        if (t == j) {
          continue;
        }
        for (int v = 0; v < n; ++v) {
          best_other[static_cast<std::size_t>(v)] = std::min(
              best_other[static_cast<std::size_t>(v)],
              dist_cols[static_cast<std::size_t>(t)][static_cast<std::size_t>(v)]);
        }
      }

      int best = medoids[static_cast<std::size_t>(j)];
      double best_cost = cost;
      std::vector<double> best_col;

      for (int c = 0; c < n_cand; ++c) {
        const int cand = members[static_cast<std::size_t>(c)];
        if (cand == medoids[static_cast<std::size_t>(j)]) {
          continue;
        }
        const std::vector<double>& col = cached_sssp(cand);
        const double c2 = trial_cost(best_other, col, n);
        if (c2 < best_cost) {
          best = cand;
          best_cost = c2;
          best_col = col;
        }
      }

      if (best != medoids[static_cast<std::size_t>(j)]) {
        medoids[static_cast<std::size_t>(j)] = best;
        dist_cols[static_cast<std::size_t>(j)].swap(best_col);
        cost = best_cost;
        assign_clusters(dist_cols, n, k, assign);
        improved = true;
      }
    }

    if (!improved) {
      break;
    }
  }

  IntegerVector out(k);
  for (int j = 0; j < k; ++j) {
    out[j] = medoids[static_cast<std::size_t>(j)];
  }
  return out;
  END_RCPP
}
