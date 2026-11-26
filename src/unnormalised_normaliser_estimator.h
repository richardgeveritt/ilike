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
  bool is_estimate() const { return true; } // (dev/plan.md 3.21 (c3))


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

  // (dev/plan.md 3.22 U0) the annealed estimate (Murray, Ghahramani & MacKay
  // 2006): each point is an annealed importance sampling run from the field
  // at theta to the auxiliary distribution at theta_hat, through `bridges`
  // distributions on the straight path between them, with `bridge_sweeps`
  // sweeps at each. 0 bridges: the importance sampling estimate.
  void set_bridges(size_t bridges_in, size_t bridge_sweeps_in);

  // (dev/plan.md 3.22 U3b) normalise = smc(particles, steps, sweeps): log
  // Z(theta_hat) estimated once, by SMC on the field from its base measure
  // (where Z = 1) along beta_t = (t/T) theta_hat, and subtracted from every
  // estimate, so that importance sampling and SMC report the evidence itself
  // (Everitt, Johansen, Rowing & Evdemon-Hogan 2017, section 2.6). The SMC
  // has its own generator, seeded from seed_value.
  void set_normalise(size_t particles_in, size_t steps_in, size_t sweeps_in, size_t seed_value);
  double get_log_normaliser_at() const { return this->log_normaliser_at; }

  // (dev/plan.md 3.22 U3b) replicates = R: R independent estimates, one per
  // replicate; the value is the sum of their logs. Under data tempering
  // (variable: the sequence's fraction), the value counts the replicates in
  // the target, and each replicate's estimate is drawn when it is first
  // needed and then kept: so between two steps of the sequence a particle's
  // value changes by the log estimates of the replicates added, drawn fresh
  // at its current parameters, which makes the SMC's incremental weight a
  // random weight (Everitt, Johansen, Rowing & Evdemon-Hogan 2017,
  // Algorithm 1) and the target at whole replicates the exact posterior.
  void set_replicates(size_t replicates_in) { this->replicates = std::max(replicates_in, size_t(1)); }
  void set_data_tempering(const std::string &variable_in) { this->tempering_variable = variable_in; }

  // log of the mean of exp(log_weights) (also used by the exchange move)
  static double log_mean_exp(const std::vector<double> &log_weights);

protected:

  friend UnnormalisedNormaliserEstimatorOutput;

  void specific_change_data(Data* new_data);

  // a draw of the field at these parameters
  arma::mat draw(const Parameters &parameters) const;
  // log of (1/M) sum_m gamma(u_m | theta_hat)/gamma(u_m | theta)
  double log_estimate(const std::vector<arma::mat> &fields, const Parameters &parameters) const;
  // one replicate's log estimate: M fresh draws at the parameters (annealed
  // with bridges), less the estimate of log Z(theta_hat) (normalise)
  double draw_replicate_log_estimate(const Parameters &parameters) const;
  // the number of replicates in the target at these parameters
  size_t active_replicates(const Parameters &parameters) const;
  // the log weight of one annealed run started from field (drawn at parameters)
  double annealed_log_weight(const arma::mat &field, const Parameters &parameters) const;

  std::string variable;
  UnnormalisedPotentialPtr potential;
  UnnormalisedExactFunctionPtr exact;   // (NULL: sweeps)
  UnnormalisedSweepFunctionPtr sweep;
  size_t sweeps;
  size_t points;
  std::vector<std::string> normaliser;  // the parameters the normaliser depends on
  std::vector<double> at;               // theta_hat, one value per parameter in normaliser
  size_t bridges = 0;
  size_t bridge_sweeps = 1;
  double log_normaliser_at = 0.0;   // (normalise: the estimate of log Z(theta_hat); 0 otherwise)
  size_t replicates = 1;
  std::string tempering_variable;   // (data tempering: the sequence's variable; empty otherwise)

  void make_copy(const UnnormalisedNormaliserEstimator &another);
};

}

#endif
