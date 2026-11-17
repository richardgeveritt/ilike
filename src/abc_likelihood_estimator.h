#ifndef ABCLIKELIHOODESTIMATOR_H
#define ABCLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <memory>

#include "likelihood_estimator.h"
#include "parameters.h"
#include "abc_kernel_factor.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

class ABCLikelihoodEstimatorOutput;
class IndependentProposalKernel;

class ABCLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  ABCLikelihoodEstimator();
  
  ABCLikelihoodEstimator(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         ABCKernelFactor* abc_kernel_in,
                         bool smcfixed_flag_in);
  
  virtual ~ABCLikelihoodEstimator();
  
  ABCLikelihoodEstimator(const ABCLikelihoodEstimator &another);
  
  void operator=(const ABCLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;
  
  LikelihoodEstimatorOutput* initialise();

  // the distance is found without random numbers (the simulation is done
  // when the point is simulated); summary: the distance (dev/plan.md 2.6a)
  bool supports_flat_generations() const { return true; }
  bool supports_parallel_chunks() const { return true; }
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);
  
private:
  
  void specific_change_data(Data* new_data);
  
  friend ABCLikelihoodEstimatorOutput;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<ABCKernelFactor> abc_kernel;
  
  std::shared_ptr<HighFive::File> h5_file;
  std::string h5_file_path;

  // values written by ABCLikelihoodEstimatorOutput, flushed per generation
  std::vector<double> pending_log_likelihoods;
  
  void make_copy(const ABCLikelihoodEstimator &another);
  
};
}

#endif
