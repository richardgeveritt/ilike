#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ilike_hdf5_utils.h"
#include <memory>

namespace ilike
{

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

  virtual ~ExactLikelihoodEstimator();
  
  ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another);
  
  void operator=(const ExactLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  LikelihoodEstimatorOutput* initialise();

  // the model is evaluated without random numbers; a summary exists when
  // smcfixed (dev/plan.md 2.6a)
  bool supports_flat_generations() const { return this->smcfixed_flag; }
  bool supports_restore_from_summary() const { return this->smcfixed_flag; }
  bool supports_reinitialise() const { return true; }
  bool supports_parallel_chunks() const { return true; }
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);

protected:
  
  void specific_change_data(Data* new_data);
  
  friend ExactLikelihoodEstimatorOutput;
  double evaluate(const Parameters &parameters) const;
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &parameters) const;
  
  double subsample_evaluate(const Parameters &parameters) const;
  arma::mat subsample_evaluate_gradient(const std::string &variable,
                                        const Parameters &parameters) const;

  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<DistributionFactor>> numerator_distribution_factors;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<LikelihoodFactor>> numerator_likelihood_factors;
  
  // not stored here
  std::vector<IndependentProposalKernel*> numerator_proposals;
  
  // not stored here
  std::vector<ProposalKernel*> numerator_likelihood_proposals;
  
  // Stored here.
  //ExactLikelihoodEstimatorOutput* output;
  
  std::shared_ptr<HighFive::File> h5_file;
  
  std::string h5_file_path;
  
  std::vector<double> pending_log_likelihoods;
  
  void make_copy(const ExactLikelihoodEstimator &another);
  
};
}

#endif
