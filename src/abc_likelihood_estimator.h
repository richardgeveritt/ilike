#ifndef ABCLIKELIHOODESTIMATOR_H
#define ABCLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "parameters.h"
#include "abc_kernel_factor.h"

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
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);

private:

  friend ABCLikelihoodEstimatorOutput;

  // Stored here.
  ABCKernelFactor* abc_kernel;
  
  std::ofstream log_likelihood_file_stream;

  void make_copy(const ABCLikelihoodEstimator &another);

};

#endif
