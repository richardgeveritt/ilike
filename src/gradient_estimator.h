#ifndef GRADIENTESTIMATOR_H
#define GRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include "particle.h"
#include "data_subsampler.h"
#include "distributions.h"

class ReinforceGradientEstimator;
class GradientEstimatorOutput;
class ProposalKernel;
class DirectGradientEstimatorOutput;

class GradientEstimator
{

public:

  GradientEstimator();
  virtual ~GradientEstimator();
  
  GradientEstimator(const ProposalKernel* proposal_in);

  GradientEstimator(const GradientEstimator &another);

  void operator=(const GradientEstimator &another);
  virtual GradientEstimator* duplicate() const=0;
  
  virtual GradientEstimatorOutput* initialise()=0;
  
  void set_proposal(const ProposalKernel* proposal_in);
  
  /*
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        Particle &particle)=0;
  
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                        Particle &particle)=0;
  
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters)=0;
  */

protected:
  
  friend ReinforceGradientEstimator;
  friend DirectGradientEstimatorOutput;
  // not stored here. A pointer to the proposal in which this gradient estimator is contained
  const ProposalKernel* proposal;
  
  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  size_t* seed;

  void make_copy(const GradientEstimator &another);

};

#endif
