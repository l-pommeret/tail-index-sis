#include <Rcpp.h>
using namespace Rcpp;

class Fenwick {
  int n;
  std::vector<int> count;
  std::vector<double> sum;
public:
  explicit Fenwick(int n_) : n(n_), count(n_ + 1, 0), sum(n_ + 1, 0.0) {}
  void add(int at, int dc, double ds) {
    for (int i = at; i <= n; i += i & -i) { count[i] += dc; sum[i] += ds; }
  }
  int prefix_count(int at) const {
    int out = 0; for (int i = at; i > 0; i -= i & -i) out += count[i]; return out;
  }
  double prefix_sum(int at) const {
    double out = 0; for (int i = at; i > 0; i -= i & -i) out += sum[i]; return out;
  }
  int kth(int target) const {
    int idx = 0, bit = 1; while ((bit << 1) <= n) bit <<= 1;
    for (; bit; bit >>= 1) {
      int next = idx + bit;
      if (next <= n && count[next] < target) { idx = next; target -= count[next]; }
    }
    return idx + 1;
  }
};

// [[Rcpp::export]]
List local_hill_cpp(NumericVector u, NumericVector y, NumericVector eval_u,
                    double h, double alpha) {
  const int n = u.size(), ng = eval_u.size();
  NumericVector estimate(ng);
  IntegerVector local_n(ng), local_k(ng);
  LogicalVector underpopulated(ng);
  std::vector<double> work;
  work.reserve(n);
  for (int g = 0; g < ng; ++g) {
    work.clear();
    for (int i = 0; i < n; ++i)
      if (std::abs(u[i] - eval_u[g]) <= h) work.push_back(std::log(y[i]));
    const int m = work.size();
    const int k = static_cast<int>(std::floor(alpha * m));
    local_n[g] = m; local_k[g] = k;
    if (alpha * m <= 1.0 || k < 1 || k >= m) {
      estimate[g] = NA_REAL; underpopulated[g] = true; continue;
    }
    std::nth_element(work.begin(), work.begin() + k, work.end(),
                     std::greater<double>());
    const double threshold = work[k];
    double total = 0.0;
    for (int i = 0; i < k; ++i) total += work[i] - threshold;
    estimate[g] = total / k;
  }
  return List::create(_["estimate"] = estimate, _["local_n"] = local_n,
                      _["local_k"] = local_k,
                      _["underpopulated"] = underpopulated);
}

// O(n log n) exact score for one coordinate, using a sliding rank window and
// Fenwick trees for local response order statistics and log-response sums.
// [[Rcpp::export]]
NumericVector score_coordinate_cpp(NumericVector x, NumericVector y, double h,
                                   double alpha, double epsilon = 0.05) {
  const int n = x.size();
  if (y.size() != n || n < 2) stop("invalid dimensions");
  std::vector<int> ox(n), oy(n), yrank(n);
  for (int i = 0; i < n; ++i) ox[i] = oy[i] = i;
  std::stable_sort(ox.begin(), ox.end(), [&](int a, int b) { return x[a] < x[b]; });
  std::stable_sort(oy.begin(), oy.end(), [&](int a, int b) { return y[a] < y[b]; });
  std::vector<double> logy_by_rank(n + 1);
  for (int r = 0; r < n; ++r) {
    yrank[oy[r]] = r + 1; logy_by_rank[r + 1] = std::log(y[oy[r]]);
  }
  Fenwick fw(n);
  int left = 0, right = -1, n_eval = 0, n_under = 0, n_grid = 0;
  double score_sum = 0.0, m_sum = 0.0, k_sum = 0.0;
  for (int pos = 0; pos < n; ++pos) {
    const double u = (pos + 1.0) / (n + 1.0);
    while (right + 1 < n && (right + 2.0) / (n + 1.0) <= u + h) {
      ++right; int id = ox[right], rr = yrank[id]; fw.add(rr, 1, logy_by_rank[rr]);
    }
    while (left <= right && (left + 1.0) / (n + 1.0) < u - h) {
      int id = ox[left], rr = yrank[id]; fw.add(rr, -1, -logy_by_rank[rr]); ++left;
    }
    if (u < epsilon || u > 1.0 - epsilon) continue;
    ++n_grid;
    const int m = right - left + 1;
    const int k = static_cast<int>(std::floor(alpha * m));
    double estimate = 0.0;
    if (alpha * m <= 1.0 || k < 1 || k >= m) {
      ++n_under;
    } else {
      const int bottom = m - k;
      const int threshold_rank = fw.kth(bottom);
      const double top_sum = fw.prefix_sum(n) - fw.prefix_sum(threshold_rank);
      const double threshold = logy_by_rank[threshold_rank];
      estimate = top_sum / k - threshold;
    }
    if (!(alpha * m <= 1.0 || k < 1 || k >= m)) {
      score_sum += estimate; ++n_eval;
    }
    m_sum += m; k_sum += k;
  }
  return NumericVector::create(_["score"] = n_eval ? score_sum / n_eval : R_PosInf,
    _["n_eval"] = n_eval, _["under_rate"] = static_cast<double>(n_under) / n_grid,
    _["mean_local_n"] = m_sum / n_grid, _["mean_local_k"] = k_sum / n_grid);
}
