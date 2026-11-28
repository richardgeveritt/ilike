#ifndef NORMALSOURCE_H
#define NORMALSOURCE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <atomic>
#include "distributions.h"

namespace ilike
{

// (dev/plan.md 3.22 U9 R3b, C306) Crank-Nicolson correlated estimates from an
// ensemble Kalman filter (Drovandi, Everitt, Golightly & Prangle 2022): the
// standard normals behind a member's draws at one step are stored (u), so
// that a proposal can move them all a little. While the filter draws member
// `column` at a step it installs a NormalSource on its thread; the built-in
// Gaussian kernels (the norm/mvnorm prior, linear_gaussian and
// nonlinear_gaussian transitions, Gaussian measurement perturbations) then
// take their standard normals from it instead of their generator. The store
// is D x N, one column per member, and gains rows (fresh N(0, 1) values, from
// `fresh`) the first time a member needs more. With no source installed (all
// other runs) the kernels draw exactly as before.
class NormalSource
{
public:
  arma::mat* store = NULL;
  arma::uword column = 0;
  arma::uword next = 0;
  RandomNumberGenerator* fresh = NULL;

  arma::colvec take(arma::uword d)
  {
    if (this->next + d > this->store->n_rows)
    {
      arma::uword old_rows = this->store->n_rows;
      this->store->resize(this->next + d, this->store->n_cols);
      for (arma::uword j = 0; j < this->store->n_cols; ++j)
        for (arma::uword r = old_rows; r < this->store->n_rows; ++r)
          (*this->store)(r, j) = rnorm(*this->fresh);
    }
    arma::colvec z = this->store->submat(this->next, this->column, this->next + d - 1, this->column);
    this->next = this->next + d;
    return z;
  }
};

// the source installed on this thread (NULL: none), and the number of
// filters, on any thread, that may install one: while it is 0 (every run
// but a correlated filter's) the kernels read only this plain counter, not
// the thread-local pointer (whose read costs a function call on macOS,
// about 1% of a particle filter's time, measured)
extern thread_local NormalSource* installed_normal_source;
extern std::atomic<int> normal_sources_in_use;

// held by a correlated filter while it draws its members or particles
struct NormalSourcesInUse
{
  NormalSourcesInUse() { normal_sources_in_use.fetch_add(1, std::memory_order_relaxed); }
  ~NormalSourcesInUse() { normal_sources_in_use.fetch_sub(1, std::memory_order_relaxed); }
  NormalSourcesInUse(const NormalSourcesInUse &) = delete;
  NormalSourcesInUse& operator=(const NormalSourcesInUse &) = delete;
};

inline NormalSource* current_normal_source()
{
  if (__builtin_expect(normal_sources_in_use.load(std::memory_order_relaxed) == 0, 1))
    return NULL;
  return installed_normal_source;
}

// As rmvnorm_using_lower_chol and rmvnorm_using_chol (distributions.h), and
// rnorm for an i.i.d. block, with the standard normals from the installed
// source when there is one.
inline arma::colvec source_rmvnorm_using_lower_chol(RandomNumberGenerator &rng,
                                                    const arma::colvec &mu,
                                                    const arma::mat &lower_chol,
                                                    double sd_scale)
{
  NormalSource* source = current_normal_source();
  if (__builtin_expect(source == NULL, 1))
    return rmvnorm_using_lower_chol(rng, mu, lower_chol, sd_scale);
  return mu + sd_scale * (lower_chol * source->take(lower_chol.n_cols));
}

inline arma::colvec source_rmvnorm_using_chol(RandomNumberGenerator &rng,
                                              const arma::colvec &mu,
                                              const arma::mat &chol)
{
  NormalSource* source = current_normal_source();
  if (__builtin_expect(source == NULL, 1))
    return rmvnorm_using_chol(rng, mu, chol);
  return mu + chol.t() * source->take(chol.n_cols);
}

inline arma::mat source_rnorm(RandomNumberGenerator &rng, std::pair<size_t, size_t> dimensions)
{
  NormalSource* source = current_normal_source();
  if (__builtin_expect(source == NULL, 1))
    return rnorm(rng, dimensions);
  return arma::reshape(source->take(dimensions.first * dimensions.second), dimensions.first, dimensions.second);
}

}

#endif
