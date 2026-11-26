#ifndef UNNORMALISEDNORMALISERESTIMATOROUTPUT_H
#define UNNORMALISEDNORMALISERESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "likelihood_estimator_output.h"

namespace ilike
{

class UnnormalisedNormaliserEstimator;

// The output of UnnormalisedNormaliserEstimator: the fields drawn at this
// particle's parameters (its auxiliary variables), and the estimate.
class UnnormalisedNormaliserEstimatorOutput : public LikelihoodEstimatorOutput
{
public:

  UnnormalisedNormaliserEstimatorOutput();
  UnnormalisedNormaliserEstimatorOutput(UnnormalisedNormaliserEstimator* estimator_in);
  virtual ~UnnormalisedNormaliserEstimatorOutput();

  UnnormalisedNormaliserEstimatorOutput(const UnnormalisedNormaliserEstimatorOutput &another);
  void operator=(const UnnormalisedNormaliserEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void simulate();
  void simulate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

  LikelihoodEstimator* get_likelihood_estimator() const;
  arma::mat get_gradient_of_log(const std::string &variable, const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable, const Parameters &x);

  // (dev/plan.md 3.21) the draws are made at the first evaluation, so that
  // a Metropolis-Hastings proposal rejected early (e.g. outside the prior's
  // range) never makes them
  bool simulation_pending() const { return this->pending; }
  void simulate_when_needed(const Parameters &parameters);

  void forget_you_were_already_written_to_file();
  void close_ofstreams();
  void write_to_file(const std::string &directory_name, const std::string &index = "");

protected:

  // Stored in ModelAndAlgorithm.
  UnnormalisedNormaliserEstimator* estimator;

  // one log estimate per replicate (dev/plan.md 3.22 U0, U3b), made from
  // fresh draws when first needed (NaN until then)
  std::vector<double> replicate_log_estimates;
  bool pending;
  void draw_if_pending(const Parameters &parameters);
  // the sum of the log estimates of the replicates in the target, drawing
  // those not drawn yet
  double current_log_estimate(const Parameters &parameters);
  double log_likelihood_smcfixed_part;

  void make_copy(const UnnormalisedNormaliserEstimatorOutput &another);
};

}

#endif
