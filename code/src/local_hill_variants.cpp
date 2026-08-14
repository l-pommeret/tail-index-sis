// Alternative functionals of the same local-Hill envelope, all obtained from
// one O(n log n) sliding-rank pass -- the pass of score_coordinate_cpp in
// local_hill.cpp, whose "mean" output is reproduced here bit for bit as a
// control.  The point of the variants is that the proposed screen ranks
// coordinates by the *level* of the envelope, which carries a large bias
// common to every coordinate; an inactive coordinate has a *flat* envelope and
// an active one a decreasing envelope, so functionals of the envelope's
// variation cancel that common level and may separate better.
//
//   mean   average local Hill over the grid            (current score)
//   min    smallest local Hill over the grid
//   range  max - min
//   sd     standard deviation over the grid
//   slope  OLS slope of the local Hill on u
//
// Kept in its own translation unit so that rebuilding it never disturbs a
// campaign already running against local_hill.cpp.
#include <Rcpp.h>
using namespace Rcpp;

namespace variants {
class Fenwick {
  int n;
  std::vector<int> count;
  std::vector<double> sum;
public:
  explicit Fenwick(int n_) : n(n_), count(n_ + 1, 0), sum(n_ + 1, 0.0) {}
  void add(int at, int dc, double ds) {
    for (int i = at; i <= n; i += i & -i) { count[i] += dc; sum[i] += ds; }
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
}

// [[Rcpp::export]]
NumericVector score_variants_cpp(NumericVector x, NumericVector y, double h,
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
  variants::Fenwick fw(n);
  int left = 0, right = -1, n_eval = 0;
  double s_e = 0.0, s_ee = 0.0, s_u = 0.0, s_uu = 0.0, s_ue = 0.0;
  double lo = R_PosInf, hi = R_NegInf;
  for (int pos = 0; pos < n; ++pos) {
    const double u = (pos + 1.0) / (n + 1.0);
    while (right + 1 < n && (right + 2.0) / (n + 1.0) <= u + h) {
      ++right; int id = ox[right], rr = yrank[id]; fw.add(rr, 1, logy_by_rank[rr]);
    }
    while (left <= right && (left + 1.0) / (n + 1.0) < u - h) {
      int id = ox[left], rr = yrank[id]; fw.add(rr, -1, -logy_by_rank[rr]); ++left;
    }
    if (u < epsilon || u > 1.0 - epsilon) continue;
    const int m = right - left + 1;
    const int k = static_cast<int>(std::floor(alpha * m));
    if (alpha * m <= 1.0 || k < 1 || k >= m) continue;
    const int threshold_rank = fw.kth(m - k);
    const double estimate = (fw.prefix_sum(n) - fw.prefix_sum(threshold_rank)) / k
                          - logy_by_rank[threshold_rank];
    ++n_eval;
    s_e += estimate; s_ee += estimate * estimate;
    s_u += u; s_uu += u * u; s_ue += u * estimate;
    if (estimate < lo) lo = estimate;
    if (estimate > hi) hi = estimate;
  }
  if (n_eval == 0)
    return NumericVector::create(_["mean"] = R_PosInf, _["min"] = R_PosInf,
      _["range"] = 0.0, _["sd"] = 0.0, _["slope"] = 0.0, _["n_eval"] = 0);
  const double mean = s_e / n_eval;
  const double var = std::max(0.0, s_ee / n_eval - mean * mean);
  const double ubar = s_u / n_eval;
  const double suu = s_uu - n_eval * ubar * ubar;
  const double sue = s_ue - n_eval * ubar * mean;
  return NumericVector::create(_["mean"] = mean, _["min"] = lo,
    _["range"] = hi - lo, _["sd"] = std::sqrt(var),
    _["slope"] = suu > 0 ? sue / suu : 0.0, _["n_eval"] = n_eval);
}
