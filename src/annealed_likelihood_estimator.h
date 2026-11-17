#ifndef ANNEALEDLIKELIHOODESTIMATOR_H
#define ANNEALEDLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <memory>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

class AnnealedLikelihoodEstimatorOutput;
class IndependentProposalKernel;

class AnnealedLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  AnnealedLikelihoodEstimator();
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              PowerFunctionPtr function_power_in,
                              const std::string &power_variable_in,
                              bool smcfixed_flag_in);
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              double constant_power_in,
                              bool smcfixed_flag_in);
  
  virtual ~AnnealedLikelihoodEstimator();
  
  AnnealedLikelihoodEstimator(const AnnealedLikelihoodEstimator &another);
  
  void operator=(const AnnealedLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // restorable when the inner estimator is (dev/plan.md 2.6a, G3)
  bool supports_restore_from_summary() const;
  bool supports_flat_generations() const;
  bool supports_parallel_chunks() const;
  bool supports_parallel_copies() const;
  void redirect_random_numbers(RandomNumberGenerator* rng_in,
                               size_t* seed_in);

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);

protected:
  
  void specific_change_data(Data* new_data);
  
  friend AnnealedLikelihoodEstimatorOutput;
  
  PowerFunctionPtr function_power;
  std::string power_variable;
  double constant_power;
  bool use_constant;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<LikelihoodEstimator> estimator;
  
  std::shared_ptr<HighFive::File> h5_file;
  std::string h5_file_path;

  // values written by AnnealedLikelihoodEstimatorOutput, flushed per generation
  std::vector<double> pending_log_likelihoods;
  
  // Stored here.
  //AnnealedLikelihoodEstimatorOutput* output;
  
  void make_copy(const AnnealedLikelihoodEstimator &another);
  
};
}

#endif
