#include "covariance_tools.h"
#include "particles.h"
#include "ensemble.h"
#include "move_output.h"

namespace ilike
{

bool weighted_covariance(const arma::mat &points,
                         const arma::colvec &normalised_weights,
                         arma::rowvec &mean,
                         arma::mat &covariance)
{
  if (points.n_rows < 2)
    return false;
  mean = normalised_weights.t()*points;
  arma::mat centred = points.each_row() - mean;
  double denominator = 1.0 - arma::accu(arma::square(normalised_weights));
  if (denominator <= 0.0)
    return false; // a single particle carries all the weight
  covariance = (centred.t()*(centred.each_col() % normalised_weights))/denominator;
  covariance = (covariance + covariance.t())/2.0;
  return true;
}

arma::mat shrunk_covariance(const arma::mat &centred,
                            const arma::colvec &normalised_weights,
                            const arma::mat &covariance,
                            double *intensity)
{
  // With S = Sb/c, Sb = sum_k w_k z_k z_k' (z_k the centred points) and
  // c = 1 - sum w_k^2, the sampling variance of entry (i, j) is estimated by
  // v_ij = sum_k w_k^2 (z_ki z_kj - Sb_ij)^2 / c^2. Summed over all (i, j):
  // sum_k w_k^2 ((z_k'z_k)^2 - 2 z_k' Sb z_k) + (sum_k w_k^2) ||Sb||^2, so no
  // p x p x n array is formed; the diagonal's share is subtracted.
  const arma::colvec &w = normalised_weights;
  arma::colvec w2 = arma::square(w);
  double c = 1.0 - arma::accu(w2);
  arma::mat biased = covariance*c;
  arma::colvec squared_norms = arma::sum(arma::square(centred), 1);
  arma::colvec quadratic = arma::sum((centred*biased) % centred, 1);
  double all = arma::accu(w2 % (arma::square(squared_norms) - 2.0*quadratic)) + arma::accu(w2)*arma::accu(arma::square(biased));
  arma::rowvec diagonal_biased = biased.diag().t();
  arma::mat squared = arma::square(centred);
  squared.each_row() -= diagonal_biased;
  arma::mat diagonal_terms = arma::square(squared);
  double diagonal = arma::accu(diagonal_terms.each_col() % w2);
  double off_variance = (all - diagonal)/(c*c);
  double off_squares = arma::accu(arma::square(covariance)) - arma::accu(arma::square(covariance.diag()));
  double lambda = 1.0;
  if (off_squares > 0.0)
    lambda = std::min(1.0, std::max(0.0, off_variance/off_squares));
  if (intensity != NULL)
    *intensity = lambda;
  arma::mat shrunk = (1.0 - lambda)*covariance;
  shrunk.diag() = covariance.diag();
  return shrunk;
}

bool conditional_moments(const arma::mat &points,
                         const arma::colvec &normalised_weights,
                         size_t block_size,
                         arma::rowvec &block_mean,
                         arma::mat &block_covariance)
{
  arma::rowvec mean;
  arma::mat covariance;
  if ((block_size == 0) || (points.n_cols < block_size) ||
      !weighted_covariance(points, normalised_weights, mean, covariance))
    return false;
  // the block, and the other columns that vary
  std::vector<arma::uword> keep;
  for (arma::uword k = 0; k < points.n_cols; ++k)
  {
    double variance = covariance(k, k);
    // (a column all points share has a variance of rounding size only)
    bool varies = std::isfinite(variance) && (variance > 1e-12*(1.0 + mean[k]*mean[k]));
    if (k < block_size)
    {
      if (!varies)
        return false;
      keep.push_back(k);
    }
    else if (varies)
      keep.push_back(k);
  }
  arma::uvec kept(keep);
  arma::mat kept_covariance = covariance(kept, kept);
  if (!kept_covariance.is_sympd())
  {
    arma::mat centred = points.cols(kept);
    centred.each_row() -= mean(kept).t();
    kept_covariance = shrunk_covariance(centred, normalised_weights, kept_covariance);
    kept_covariance = (kept_covariance + kept_covariance.t())/2.0;
    if (!kept_covariance.is_sympd())
      return false;
  }
  block_mean = mean.head(block_size);
  if (kept.n_elem == block_size)
  {
    block_covariance = kept_covariance;
    return true;
  }
  // ((Sigma^-1)_xx)^-1
  arma::mat precision;
  if (!arma::inv_sympd(precision, kept_covariance))
    return false;
  arma::mat block_precision = precision.submat(0, 0, block_size - 1, block_size - 1);
  block_precision = (block_precision + block_precision.t())/2.0;
  if (!arma::inv_sympd(block_covariance, block_precision))
    return false;
  block_covariance = (block_covariance + block_covariance.t())/2.0;
  return true;
}

bool conditional_variances(const arma::mat &points,
                           const arma::colvec &normalised_weights,
                           arma::colvec &variances)
{
  arma::rowvec mean;
  arma::mat covariance;
  if (!weighted_covariance(points, normalised_weights, mean, covariance))
    return false;
  std::vector<arma::uword> keep;
  for (arma::uword k = 0; k < points.n_cols; ++k)
    if (std::isfinite(covariance(k, k)) && (covariance(k, k) > 1e-12*(1.0 + mean[k]*mean[k])))
      keep.push_back(k);
  variances = arma::colvec(points.n_cols, arma::fill::value(arma::datum::nan));
  if (keep.empty())
    return true;
  arma::uvec kept(keep);
  arma::mat kept_covariance = covariance(kept, kept);
  if (!kept_covariance.is_sympd())
  {
    arma::mat centred = points.cols(kept);
    centred.each_row() -= mean(kept).t();
    kept_covariance = shrunk_covariance(centred, normalised_weights, kept_covariance);
    kept_covariance = (kept_covariance + kept_covariance.t())/2.0;
    if (!kept_covariance.is_sympd())
      return false;
  }
  arma::mat precision;
  if (!arma::inv_sympd(precision, kept_covariance))
    return false;
  for (arma::uword k = 0; k < kept.n_elem; ++k)
    variances[kept[k]] = 1.0/precision(k, k);
  return true;
}

// the variables' names (the given ones first) and offsets, from one point
static void layout(const Parameters &first,
                   std::vector<std::string> &names,
                   std::vector<size_t> &offsets)
{
  std::vector<std::string> all = first.get_vector_variables();
  for (const std::string &name : all)
    if (std::find(names.begin(), names.end(), name) == names.end())
      names.push_back(name);
  offsets.clear();
  size_t offset = 0;
  for (const std::string &name : names)
  {
    offsets.push_back(offset);
    offset = offset + first[name].n_elem;
  }
}

arma::mat population_points(const Particles &generation,
                            std::vector<std::string> &names,
                            std::vector<size_t> &offsets)
{
  if (generation.is_flat())
  {
    if (generation.flat_points.size() == 0)
      return arma::mat();
    layout(generation.flat_points.parameters_at(0), names, offsets);
    return generation.flat_points.matrix_of(names);
  }
  size_t n = generation.particles.size();
  if (n == 0)
    return arma::mat();
  layout(generation.particles[0]->back().parameters, names, offsets);
  arma::mat points;
  for (size_t i = 0; i < n; ++i)
  {
    arma::rowvec point = generation.particles[i]->back().parameters.get_rowvec(names);
    if (i == 0)
      points.set_size(n, point.n_elem);
    points.row(i) = point;
  }
  return points;
}

arma::mat ensemble_points(const Ensemble &ensemble,
                          std::vector<std::string> &names,
                          std::vector<size_t> &offsets)
{
  size_t n = ensemble.members.size();
  if (n == 0)
    return arma::mat();
  layout(ensemble.members[0]->back().parameters, names, offsets);
  arma::mat points;
  for (size_t i = 0; i < n; ++i)
  {
    arma::rowvec point = ensemble.members[i]->back().parameters.get_rowvec(names);
    if (i == 0)
      points.set_size(n, point.n_elem);
    points.row(i) = point;
  }
  return points;
}

bool population_weights(const Particles &generation,
                        size_t number_of_points,
                        arma::colvec &weights)
{
  // the weights for the current target, before resampling
  const arma::colvec &log_weights = generation.unnormalised_log_weights;
  if (log_weights.n_elem == number_of_points)
  {
    weights = arma::exp(log_weights - log_weights.max());
    if (!weights.is_finite() || (arma::accu(weights) <= 0.0))
      return false;
    weights = weights/arma::accu(weights);
  }
  else
    weights = arma::colvec(number_of_points, arma::fill::value(1.0/double(number_of_points)));
  return true;
}

}

// For tests: the numerics above on given points (dev/plan.md, "conditional
// scaling of partial moves").
// [[Rcpp::export]]
Rcpp::List covariance_tools_check(const arma::mat &points, const arma::colvec &weights, int block_size)
{
  arma::colvec w = weights/arma::accu(weights);
  arma::rowvec mean;
  arma::mat covariance;
  ilike::weighted_covariance(points, w, mean, covariance);
  arma::mat centred = points.each_row() - mean;
  double lambda = 0.0;
  arma::mat shrunk = ilike::shrunk_covariance(centred, w, covariance, &lambda);
  arma::rowvec block_mean;
  arma::mat block_covariance;
  bool ok = ilike::conditional_moments(points, w, size_t(block_size), block_mean, block_covariance);
  arma::colvec variances;
  bool ok_variances = ilike::conditional_variances(points, w, variances);
  return Rcpp::List::create(Rcpp::Named("covariance") = covariance, Rcpp::Named("shrunk") = shrunk,
                            Rcpp::Named("lambda") = lambda, Rcpp::Named("ok") = ok,
                            Rcpp::Named("block_mean") = block_mean, Rcpp::Named("block_covariance") = block_covariance,
                            Rcpp::Named("ok_variances") = ok_variances, Rcpp::Named("variances") = variances);
}
