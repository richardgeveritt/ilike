#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"

class ExactLikelihoodEstimatorOutput;
class IndependentProposalKernel;
class DistributionFactor;
class LikelihoodFactor;
class ProposalKernel;

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
                           DistributionFactor* dist_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           LikelihoodFactor* llhd_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                           const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           const std::vector<DistributionFactor*> &numerator_distribution_factors_in,
                           const std::vector<LikelihoodFactor*> &numerator_likelihood_factors_in,
                           const std::vector<ProposalKernel*> &numerator_likelihood_proposals_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           IndependentProposalKernel* dist_in,
                           bool smcfixed_flag_in);
  
  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           ProposalKernel* dist_in,
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
                           Data* data_in,
                           std::vector<DistributionFactor*> numerator_distribution_factors,
                           bool smcfixed_flag_in);
  */
  
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
  
  void setup();
  void setup(const Parameters &parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
protected:
  
  void specific_change_data(Data* new_data);
  
  friend ExactLikelihoodEstimatorOutput;
  double evaluate(const Parameters &parameters) const;
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &parameters) const;
  
  double subsample_evaluate(const Parameters &parameters) const;
  arma::mat subsample_evaluate_gradient(const std::string &variable,
                                        const Parameters &parameters) const;
  //double evaluate(const Parameters &parameters);
  
  /*
  std::vector<EvaluateLogLikelihoodPtr> numerator_llhds;
  std::vector<EvaluateLogDistributionPtr> numerator_distributions;
  std::vector<EvaluateLogLikelihoodPtr> denominator_llhds;
  std::vector<EvaluateLogDistributionPtr> denominator_distributions;
  
  std::vector<EvaluateGradientLogLikelihoodPtr> gradient_numerator_llhds;
  std::vector<EvaluateGradientLogDistributionPtr> gradient_numerator_distributions;
  std::vector<EvaluateGradientLogLikelihoodPtr> gradient_denominator_llhds;
  std::vector<EvaluateGradientLogDistributionPtr> gradient_denominator_distributions;
  */
  
  // stored here
  std::vector<DistributionFactor*> numerator_distribution_factors;
  
  // stored here
  std::vector<LikelihoodFactor*> numerator_likelihood_factors;
  
  // not stored here
  std::vector<IndependentProposalKernel*> numerator_proposals;
  
  // not stored here
  std::vector<ProposalKernel*> numerator_likelihood_proposals;
  
  // Stored here.
  //ExactLikelihoodEstimatorOutput* output;
  
  std::ofstream log_likelihood_file_stream;
  
  void make_copy(const ExactLikelihoodEstimator &another);
  
};

#endif
