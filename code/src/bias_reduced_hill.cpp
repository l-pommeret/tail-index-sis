// Second-order bias-reduced local Hill.
//
// Diagnosis this addresses: with k/m near 8% the local Hill does not isolate
// the upper envelope of the conditional tail-index surface along a fibre, it
// estimates a blend.  Measured on M2 at n=5000, p=2000, every coordinate scores
// about 0.15 while the population envelope is 0.313 for an active coordinate
// and 0.500 for an inactive one: the separation is compressed from 0.186 to
// 0.008.  Shrinking k removes the blend but explodes the variance (measured:
// separation-to-noise falls from 3.50 at a=0.40 to 2.27 at a=0.50).  Bias
// reduction is the only move that removes the blend while KEEPING k large.
//
// Construction (Feuerverger-Hall / Beirlant et al.).  For the k largest
// responses in a rank window, the scaled log-spacings
//     Z_i = i * (log Y_(i) - log Y_(i+1)),   i = 1..k
// are approximately independent exponentials with
//     E[Z_i] ~ gamma + b * (i/(k+1))^(-rho),   rho <= 0 the second-order index.
// Regressing Z_i on (i/(k+1))^(-rho) by least squares gives an intercept that
// estimates gamma free of the first-order bias term; the plain Hill estimator
// is exactly the mean of the Z_i, and is returned alongside as the control, so
// the comparison isolates the bias correction and not the evaluation grid.
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector score_bias_reduced_cpp(NumericVector x, NumericVector y,
                                     double h, double alpha, double rho = -1.0,
                                     double epsilon = 0.05, int n_grid = 25) {
  const int n = x.size();
  if (y.size() != n || n < 2) stop("invalid dimensions");
  if (rho > 0) stop("rho must be non-positive");
  // rank-transform x once: the window is defined on empirical ranks
  std::vector<int> ox(n);
  for (int i = 0; i < n; ++i) ox[i] = i;
  std::stable_sort(ox.begin(), ox.end(), [&](int a, int b) { return x[a] < x[b]; });
  std::vector<double> u(n), ly(n);
  for (int r = 0; r < n; ++r) { u[ox[r]] = (r + 1.0) / (n + 1.0); }
  for (int i = 0; i < n; ++i) ly[i] = std::log(y[i]);

  double sum_hill = 0.0, sum_br = 0.0;
  int n_eval = 0;
  std::vector<double> work;
  work.reserve(n);
  for (int g = 0; g < n_grid; ++g) {
    const double center = epsilon +
      (1.0 - 2.0 * epsilon) * (g + 0.5) / n_grid;
    work.clear();
    for (int i = 0; i < n; ++i)
      if (std::abs(u[i] - center) <= h) work.push_back(ly[i]);
    const int m = work.size();
    const int k = static_cast<int>(std::floor(alpha * m));
    if (alpha * m <= 1.0 || k < 3 || k + 1 >= m) continue;
    // k+1 largest, in descending order
    std::nth_element(work.begin(), work.begin() + (k + 1), work.end(),
                     std::greater<double>());
    std::sort(work.begin(), work.begin() + (k + 1), std::greater<double>());
    // scaled log-spacings and their regression on (i/(k+1))^(-rho)
    double sz = 0.0, sx = 0.0, sxx = 0.0, sxz = 0.0;
    for (int i = 1; i <= k; ++i) {
      const double z = i * (work[i - 1] - work[i]);
      const double xi = std::pow(static_cast<double>(i) / (k + 1.0), -rho);
      sz += z; sx += xi; sxx += xi * xi; sxz += xi * z;
    }
    const double zbar = sz / k, xbar = sx / k;
    const double sxx_c = sxx - k * xbar * xbar;
    const double sxz_c = sxz - k * xbar * zbar;
    const double slope = sxx_c > 0 ? sxz_c / sxx_c : 0.0;
    sum_hill += zbar;                    // mean of the Z_i is the Hill estimator
    sum_br += zbar - slope * xbar;       // intercept: bias-reduced
    ++n_eval;
  }
  if (n_eval == 0)
    return NumericVector::create(_["hill"] = R_PosInf, _["br"] = R_PosInf,
                                 _["n_eval"] = 0);
  return NumericVector::create(_["hill"] = sum_hill / n_eval,
                               _["br"] = sum_br / n_eval,
                               _["n_eval"] = n_eval);
}
