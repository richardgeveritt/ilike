#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"

class ExactLikelihoodEstimatorOutput;
class IndependentProposalKernel;

class ExactLikelihoodEstimator : public LikelihoodEstimator
{

public:

  ExactLikelihoodEstimator();

  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           EvaluateLogLikelihoodPtr llhd_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           EvaluateLogDistributionPtr dist_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           IndependentProposalKernel* dist_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           EvaluateLogDistributionPtr prior_in,
                           EvaluateLogLikelihoodPtr llhd_in,
                           bool smcfixed_flag_in);
  
  /*
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           const Data* data_in,
                           EvaluateLogDistributionPtr prior_in,
                           EvaluateLogLikelihoodPtr llhd_in,
                           EvaluateLogDistributionPtr proposal_in,
                           bool smcfixed_flag_in);
   */

  virtual ~ExactLikelihoodEstimator();

  ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another);

  void operator=(const ExactLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

protected:

  friend ExactLikelihoodEstimatorOutput;
  double evaluate(const Parameters &parameters);
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &parameters);
  
  double subsample_evaluate(const Parameters &parameters);
  arma::mat subsample_evaluate_gradient(const std::string &variable,
                                        const Parameters &parameters);
  //double evaluate(const Parameters &parameters);
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;

  std::vector<EvaluateLogLikelihoodPtr> numerator_llhds;
  std::vector<EvaluateLogDistributionPtr> numerator_distributions;
  std::vector<EvaluateLogLikelihoodPtr> denominator_llhds;
  std::vector<EvaluateLogDistributionPtr> denominator_distributions;
  
  std::vector<EvaluateGradientLogLikelihoodPtr> gradient_numerator_llhds;
  std::vector<EvaluateGradientLogDistributionPtr> gradient_numerator_distributions;
  std::vector<EvaluateGradientLogLikelihoodPtr> gradient_denominator_llhds;
  std::vector<EvaluateGradientLogDistributionPtr> gradient_denominator_distributions;
  
  // not stored here
  std::vector<IndependentProposalKernel*> numerator_proposals;

  // Stored here.
  //ExactLikelihoodEstimatorOutput* output;

  void make_copy(const ExactLikelihoodEstimator &another);

};

#endif
