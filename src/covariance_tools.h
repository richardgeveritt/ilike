#ifndef COVARIANCETOOLS_H
#define COVARIANCETOOLS_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>

namespace ilike
{

class Particles;
class Ensemble;

// Covariance estimates for adapting proposals that move part of the
// parameters (dev/plan.md, "conditional scaling of partial moves"): a move
// that changes a block x of the parameters, the rest r fixed, explores the
// conditional distribution of x given r, whose covariance (for a Gaussian
// approximation) is Sigma_x|r = ((Sigma^-1)_xx)^-1, not Sigma_xx.

// The weighted covariance of the points (one row per point, normalised
// weights; the unbiased form, divided by 1 - sum of squared weights), and
// their mean. FALSE if a single point carries all the weight.
bool weighted_covariance(const arma::mat &points,
                         const arma::colvec &normalised_weights,
                         arma::rowvec &mean,
                         arma::mat &covariance);

// Shrinkage of a weighted covariance towards its diagonal (Ledoit & Wolf
// 2004; the intensity of Schafer & Strimmer 2005, target "D", with weights):
// (1 - lambda) S + lambda diag(S), lambda estimated as the summed sampling
// variance of the off-diagonal entries over their summed squares. centred:
// the points minus their weighted mean. With positive variances it is
// positive definite whenever lambda > 0.
arma::mat shrunk_covariance(const arma::mat &centred,
                            const arma::colvec &normalised_weights,
                            const arma::mat &covariance,
                            double *intensity = NULL);

// The mean and conditional covariance of the first block_size columns given
// the others, from weighted points: the other columns that do not vary
// (e.g. a variable all points share) are left out; the covariance is shrunk
// (shrunk_covariance) when it is not positive definite. FALSE if it cannot
// be estimated (a single point has all the weight, a column of the block does
// not vary, or the covariance is still singular): the caller keeps its
// proposal.
bool conditional_moments(const arma::mat &points,
                         const arma::colvec &normalised_weights,
                         size_t block_size,
                         arma::rowvec &block_mean,
                         arma::mat &block_covariance);

// The conditional variance of each column given all the others, from one
// precision matrix (1/(Sigma^-1)_jj): for many single-element moves at once.
// Columns that do not vary get NaN. FALSE if it cannot be estimated.
bool conditional_variances(const arma::mat &points,
                           const arma::colvec &normalised_weights,
                           arma::colvec &variances);

// All the vector variables of the points of an SMC generation (one row per
// particle), in the order of `names` (filled from the first particle when
// empty, the given ones first), and where each variable starts.
arma::mat population_points(const Particles &generation,
                            std::vector<std::string> &names,
                            std::vector<size_t> &offsets);
arma::mat ensemble_points(const Ensemble &ensemble,
                          std::vector<std::string> &names,
                          std::vector<size_t> &offsets);

// The normalised weights of an SMC generation for its current target (before
// resampling), or equal weights; FALSE if they are not usable.
bool population_weights(const Particles &generation,
                        size_t number_of_points,
                        arma::colvec &weights);

}

#endif
