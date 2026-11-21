#ifndef UNNORMALISEDNORMALISERESTIMATOR_H
#define UNNORMALISEDNORMALISERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "likelihood_estimator.h"

namespace ilike
{

// An estimate of 1/Z(theta), the normalising constant of an unnormalised
// distribution gamma(x | theta)/Z(theta) (dev/plan.md 3.20, S4): the single
// (and multiple) auxiliary variable methods of Moller et al. (2006) and
// Murray et al. (2006). With u_1, ..., u_M drawn from the distribution at
// theta, (1/M) sum_m gamma(u_m | theta_hat)/gamma(u_m | theta) has
// expectation Z(theta_hat)/Z(theta): an unbiased estimate of 1/Z(theta) up to
// the constant Z(theta_hat), which cancels in every ratio. The draws are this
// factor's auxiliary variables (pseudo-marginal: they are kept with the
// state, and their density cancels with their simulation).
typedef double (*UnnormalisedPotentialPtr)(const arma::mat &field, const Parameters &parameters, const Data &data);
typedef void (*UnnormalisedSweepFunctionPtr)(RandomNumberGenerator &rng, Parameters &parameters, const Data &data);
typedef Parameters (*UnnormalisedExactFunctionPtr)(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data);

class UnnormalisedNormaliserEstimatorOutput;

class UnnormalisedNormaliserEstimator : public LikelihoodEstimator
{
public:

  UnnormalisedNormaliserEstimator();
  UnnormalisedNormaliserEstimator(RandomNumberGenerator* rng_in,
                                  size_t* seed_in,
                                  Data* data_in,
                                  const std::string &variable_in,
                                  UnnormalisedPotentialPtr potential_in,
                                  UnnormalisedExactFunctionPtr exact_in,
                                  UnnormalisedSweepFunctionPtr sweep_in,
                                  size_t sweeps_in,
                                  size_t points_in,
                                  const std::vector<std::string> &normaliser_in,
                                  const std::vector<double> &at_in,
                                  bool smcfixed_flag_in);
  virtual ~UnnormalisedNormaliserEstimator();

  UnnormalisedNormaliserEstimator(const UnnormalisedNormaliserEstimator &another);
  void operator=(const UnnormalisedNormaliserEstimator &another);
  LikelihoodEstimator* duplicate() const;
  bool defers_simulation() const { return true; }

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  void setup();
  void setup(const Parameters &parameters);

protected:

  friend UnnormalisedNormaliserEstimatorOutput;

  void specific_change_data(Data* new_data);

  // a draw of the field at these parameters
  arma::mat draw(const Parameters &parameters) const;
  // log of (1/M) sum_m gamma(u_m | theta_hat)/gamma(u_m | theta)
  double log_estimate(const std::vector<arma::mat> &fields, const Parameters &parameters) const;

  std::string variable;
  UnnormalisedPotentialPtr potential;
  UnnormalisedExactFunctionPtr exact;   // (NULL: sweeps)
  UnnormalisedSweepFunctionPtr sweep;
  size_t sweeps;
  size_t points;
  std::vector<std::string> normaliser;  // the parameters the normaliser depends on
  std::vector<double> at;               // theta_hat, one value per parameter in normaliser

  void make_copy(const UnnormalisedNormaliserEstimator &another);
};

}

#endif
