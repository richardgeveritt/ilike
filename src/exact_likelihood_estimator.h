#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file exact_likelihood_estimator.h
   * @brief Defines the ExactLikelihoodEstimatorOutput class.
   *
   * Stores the output of a exact likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class ExactLikelihoodEstimatorOutput
   * @brief The exact likelihood estimator output class.
   */


class ExactLikelihoodEstimatorOutput;
class IndependentProposalKernel;
class DistributionFactor;
class LikelihoodFactor;
class ProposalKernel;

class ExactLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the exactlikelihoodestimator operation.
   */
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
  
  /**
   * @brief Performs the exactlikelihoodestimator operation.
   *
   * @param another The ExactLikelihoodEstimatorOutput instance to copy from.
   */
  ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another);
  
  /**
   * @brief Assignment operator for ExactLikelihoodEstimatorOutput.
   *
   * @param another The ExactLikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const ExactLikelihoodEstimator &another);
  /**
   * @brief Creates a deep copy of this ExactLikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend ExactLikelihoodEstimatorOutput;
  /**
   * @brief Evaluates.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double evaluate(const Parameters &parameters) const;
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &parameters) const;
  
  /**
   * @brief Performs the subsample evaluate operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
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
  /** @brief The numerator distribution factors. */
  std::vector<DistributionFactor*> numerator_distribution_factors;
  
  // stored here
  /** @brief The numerator likelihood factors. */
  std::vector<LikelihoodFactor*> numerator_likelihood_factors;
  
  // not stored here
  /** @brief The numerator proposals. */
  std::vector<IndependentProposalKernel*> numerator_proposals;
  
  // not stored here
  /** @brief The numerator likelihood proposals. */
  std::vector<ProposalKernel*> numerator_likelihood_proposals;
  
  // Stored here.
  //ExactLikelihoodEstimatorOutput* output;
  
  /** @brief The log likelihood file stream. */
  std::ofstream log_likelihood_file_stream;
  
  /**
   * @brief Copies the state of another ExactLikelihoodEstimatorOutput into this object.
   *
   * @param another The ExactLikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const ExactLikelihoodEstimator &another);
  
};
}

#endif
