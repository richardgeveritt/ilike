#include <algorithm>
#include <numeric>
#include "hilbert_sort.h"

namespace ilike
{

uint64_t hilbert_index(uint64_t* x, size_t n, unsigned bits)
{
  if ((n == 0) || (bits == 0))
    return 0;
  const uint64_t M = uint64_t(1) << (bits - 1);
  // Skilling's AxesToTranspose: inverse undo excess work
  for (uint64_t Q = M; Q > 1; Q >>= 1)
  {
    const uint64_t P = Q - 1;
    for (size_t i = 0; i < n; ++i)
    {
      if (x[i] & Q)
        x[0] ^= P;
      else
      {
        uint64_t t = (x[0] ^ x[i]) & P;
        x[0] ^= t;
        x[i] ^= t;
      }
    }
  }
  // Gray encode
  for (size_t i = 1; i < n; ++i)
    x[i] ^= x[i - 1];
  uint64_t t = 0;
  for (uint64_t Q = M; Q > 1; Q >>= 1)
    if (x[n - 1] & Q)
      t ^= Q - 1;
  for (size_t i = 0; i < n; ++i)
    x[i] ^= t;
  // the transposed index, most significant bits first
  uint64_t h = 0;
  for (int b = int(bits) - 1; b >= 0; --b)
    for (size_t i = 0; i < n; ++i)
      h = (h << 1) | ((x[i] >> b) & uint64_t(1));
  return h;
}

std::vector<size_t> hilbert_order(const arma::mat &points)
{
  const size_t n = points.n_rows, k = points.n_cols;
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), size_t(0));
  if ((n < 2) || (k == 0))
    return order;
  if (k == 1)
  {
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return points(a, 0) < points(b, 0); });
    return order;
  }
  if (k > 64)
    Rcpp::stop("correlated = crank_nicolson(…) in a particle filter: the Hilbert sort of states with more than 64 coordinates is not available.");
  // (12 bits: 4096 cells per coordinate, about 0.002 standard deviations
  // wide near the centre, so that particles rarely share a cell; more only
  // costs time, the index's work growing with the bits)
  const unsigned bits = unsigned(std::min<size_t>(64 / k, 12));
  const double cells = std::ldexp(1.0, int(bits));
  std::vector<uint64_t> keys(n);
  arma::rowvec mean = arma::mean(points, 0);
  arma::rowvec sd = arma::stddev(points, 0, 0);
  for (size_t j = 0; j < k; ++j)
    sd(j) = (sd(j) > 0.0) ? sd(j) : 1.0;
  uint64_t cell[64];
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < k; ++j)
    {
      double u = 1.0 / (1.0 + std::exp(-(points(i, j) - mean(j)) / sd(j)));
      double c = std::floor(u * cells);
      cell[j] = uint64_t(std::min(std::max(c, 0.0), cells - 1.0));
    }
    keys[i] = hilbert_index(cell, k, bits);
  }
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return keys[a] < keys[b]; });
  return order;
}

}
