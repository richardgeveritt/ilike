#ifndef UNNORMALISEDNORMALISERESTIMATOR_H
#define UNNORMALISEDNORMALISERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "likelihood_estimator.h"
#include "annealed_path.h"
#include "smc_path.h"

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
  // (dev/plan.md 3.22 U3c) resampling = ess(p): the M annealed runs moved
  // together, resampled when their effective sample size falls below p M
  // (as SMC's resampling = ess(p)); 0: independent runs
  void set_resampling(double threshold_in) { this->resampling_threshold = threshold_in; }

  // (dev/plan.md 3.22 U3b) normalise = smc(particles, steps, sweeps): log
  // Z(theta_hat) estimated once, by SMC on the field from its base measure
  // (where Z = 1) along beta_t = (t/T) theta_hat, and subtracted from every
  // estimate, so that importance sampling and SMC report the evidence itself
  // (Everitt, Johansen, Rowing & Evdemon-Hogan 2017, section 2.6). The SMC
  // has its own generator, seeded from seed_value.
  void set_normalise(size_t particles_in, size_t steps_in, size_t sweeps_in, size_t seed_value, double resampling_threshold_in = 1.0);
  double get_log_normaliser_at() const { return this->log_normaliser_at; }

  // (dev/plan.md 3.22 U9 R1) normaliser = estimate(smc(…)): the value is
  // -log Zhat(theta), Zhat an SMC estimate from the base measure at the
  // current parameters (the same SMC as normalise's), kept with the state as
  // any estimated factor. exp(-log Zhat) is not unbiased for 1/Z, so the
  // posterior is approximate (the noisy plug-in; the debiasing estimators of
  // U9 R2-R4 remove the bias).
  void set_plugin(size_t particles_in, size_t steps_in, size_t sweeps_in, double resampling_threshold_in);
  // (dev/plan.md 3.22 U9 R2) debias = roulette(continue = q, scale = k): each
  // replicate's value an unbiased, possibly negative, estimate of 1/Z(theta),
  // Lyne, Girolami, Atchadé, Strathmann & Simpson's (2015) geometric series
  //   1/Z = c sum_{n >= 0} (1 - cZ)^n,  c = 1/(k Zhat_0(theta)),
  // Zhat_0 an independent pilot estimate (so that cZ is near 1/k wherever theta
  // is), each (1 - cZ)^n estimated by prod_{i <= n} (1 - c Zhat_i), and the
  // series truncated at random: continued with probability q, the n-th term
  // divided by q^n. Its log absolute value is the factor's value, and its sign
  // the output's (signed pseudo-marginal MCMC).
  void set_roulette(double continue_in, double scale_in);
  // (dev/plan.md 3.22 U9 R3a) debias = block_poisson(blocks = lambda, mean =
  // m, lower = a, correlated): the block-Poisson estimate of exp(-nu Z(theta))
  // (Yang, Quiroz, Kohn & Sisson, arXiv 2210.02734, Definition 1), with nu an
  // auxiliary variable drawn from Exp(Zhat_P(theta)): the output holds the
  // block seeds and nu, and the value is
  //   log |Lhat_B| - log Zhat_P + nu Zhat_P,
  // so that a Metropolis(-Hastings) move on theta is the paper's Algorithm 1
  // (UnnormalisedNormaliserEstimatorOutput::block_poisson_log_estimate).
  void set_block_poisson(size_t blocks_in, double mean_in, double lower_in, bool correlated_in);
  bool may_be_negative() const { return this->roulette || this->block_poisson; }

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
  // (dev/plan.md 3.22 U4) the general engine, for a distribution written as
  // a block: its draws (draw_sweeps iterations of its kernel, from a base
  // draw or the latent field), its annealed paths (one particle, and the
  // population's), and the SMC of normalise. Owned; any may be NULL.
  void set_engine(SMCPath* single_in, SMCPath* population_in, SMCPath* normalise_in, size_t draw_sweeps_in);
  void set_data_tempering(const std::string &variable_in) { this->tempering_variable = variable_in; }

  // log of the mean of exp(log_weights) (also used by the exchange move)
  static double log_mean_exp(const std::vector<double> &log_weights);

protected:

  friend UnnormalisedNormaliserEstimatorOutput;

  void specific_change_data(Data* new_data);

  // log Zhat(end) relative to the base measure: an SMC from the base measure
  // (the general engine's, or the fast engine's along beta from 0 to end)
  double base_log_evidence(RandomNumberGenerator &rng, const Parameters &end, size_t particles, size_t steps,
                           size_t sweeps_in, double resampling_threshold_in) const;
  // a draw of the field at these parameters
  arma::mat draw(const Parameters &parameters) const;
  // one replicate's log estimate: M fresh draws at the parameters (annealed
  // with bridges), less the estimate of log Z(theta_hat) (normalise)
  double draw_replicate_log_estimate(const Parameters &parameters, double *sign = NULL) const;
  // the number of replicates in the target at these parameters
  size_t active_replicates(const Parameters &parameters) const;
  // the log weight of one annealed run started from field (drawn at parameters)
  double annealed_log_weight(const arma::mat &field, const Parameters &parameters) const;
  // (dev/plan.md 3.22 U3c) the path to theta_hat, through the shared fast engine
  FieldPath path() const;
  Parameters at_parameters() const;

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
  double resampling_threshold = 0.0;
  double log_normaliser_at = 0.0;   // (normalise: the estimate of log Z(theta_hat); 0 otherwise)
  size_t replicates = 1;
  std::string tempering_variable;   // (data tempering: the sequence's variable; empty otherwise)
  ClonePtr<SMCPath> engine_single, engine_population, engine_normalise;
  bool plugin = false;
  size_t plugin_particles = 200, plugin_steps = 100, plugin_sweeps = 1;
  double plugin_resampling = 1.0;
  bool roulette = false;
  double roulette_continue = 0.5, roulette_scale = 1.0;
  bool block_poisson = false;
  size_t block_poisson_blocks = 10;
  double block_poisson_mean = 1.0, block_poisson_lower = -11.0;
  bool block_poisson_correlated = true;
  // (block-Poisson) one SMC estimate of log Z(theta) relative to the base
  // measure, from the generator given (the plug-in's SMC)
  double plugin_log_z(RandomNumberGenerator &rng, const Parameters &parameters) const;

  void make_copy(const UnnormalisedNormaliserEstimator &another);
};

}

#endif
