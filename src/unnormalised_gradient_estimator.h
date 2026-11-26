#ifndef UNNORMALISEDGRADIENTESTIMATOR_H
#define UNNORMALISEDGRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "likelihood_estimator.h"
#include "unnormalised_normaliser_estimator.h"

namespace ilike
{

class UnnormalisedGradientEstimatorOutput;

// The gradient of log 1/Z(beta) for an unnormalised distribution
// gamma(x | beta)/Z(beta) (dev/plan.md 3.22 U1b; Alquier, Friel, Everitt &
// Boland 2016, eq. (7)): grad log Z(beta) = E_beta[d log gamma(u | beta)/d beta],
// estimated by the mean over N fields u_i drawn at beta. A factor of its own
// (x_log_normaliser, from gradient = estimate(...)) whose *value* is 0: the
// normalising constant is not known, and the moves that may change beta
// handle it in their acceptance ratios (exchange), or never evaluate factors
// (unadjusted moves); only its gradient is used, by gradient-based proposals.
// The fields are drawn at the first gradient asked for, and kept with the
// particle (so that a Metropolis-Hastings proposal's density at the current
// point uses the same fields as the proposal did: MALA-exchange).
class UnnormalisedGradientEstimator : public LikelihoodEstimator
{
public:

  bool is_estimate() const { return true; }

  UnnormalisedGradientEstimator();

  UnnormalisedGradientEstimator(RandomNumberGenerator* rng_in,
                                size_t* seed_in,
                                Data* data_in,
                                const std::string &variable_in,
                                UnnormalisedPotentialPtr statistic_in,
                                UnnormalisedExactFunctionPtr exact_in,
                                UnnormalisedSweepFunctionPtr sweep_in,
                                size_t sweeps_in,
                                size_t points_in,
                                const std::string &parameter_in);

  virtual ~UnnormalisedGradientEstimator();

  UnnormalisedGradientEstimator(const UnnormalisedGradientEstimator &another);
  void operator=(const UnnormalisedGradientEstimator &another);
  LikelihoodEstimator* duplicate() const;

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  void setup();
  void setup(const Parameters &parameters);

protected:

  friend UnnormalisedGradientEstimatorOutput;

  void specific_change_data(Data* new_data);

  // -(1/N) sum_i d log gamma(u_i | beta)/d beta, u_i drawn at the parameters
  double estimate(const Parameters &parameters) const;

  std::string variable;                  // the field
  UnnormalisedPotentialPtr statistic;    // d log gamma(field | beta)/d beta
  UnnormalisedExactFunctionPtr exact;    // (NULL: sweeps)
  UnnormalisedSweepFunctionPtr sweep;
  size_t sweeps;
  size_t points;
  std::string parameter;                 // beta
  size_t replicates = 1;                 // (replicates = R: the gradient of R log 1/Z)
public:
  void set_replicates(size_t replicates_in) { this->replicates = std::max(replicates_in, size_t(1)); }
protected:

  void make_copy(const UnnormalisedGradientEstimator &another);
};

}

#endif
