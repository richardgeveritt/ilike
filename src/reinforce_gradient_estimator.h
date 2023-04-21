#ifndef REINFORCEGRADIENTESTIMATOR_H
#define REINFORCEGRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator.h"
#include "gaussian_random_walk_proposal_kernel.h"

class GradientEstimatorOutput;

class ReinforceGradientEstimator : public GradientEstimator
{

public:

  ReinforceGradientEstimator();
  
  ReinforceGradientEstimator(DataSubsampler* subsampler_in);

  virtual ~ReinforceGradientEstimator();

  ReinforceGradientEstimator(const ReinforceGradientEstimator &another);

  void operator=(const ReinforceGradientEstimator &another);
  GradientEstimator* duplicate() const;
  
  GradientEstimatorOutput* initialise() const;

  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle);
  
  /*
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters);
  */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle);
  
  /*
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters);
  */
  
protected:

  void make_copy(const ReinforceGradientEstimator &another);
  
  // can be generalised, but let's not worry about that for now
  GaussianRandomWalkProposalKernel gaussian_proposal;
  
  size_t num_points;
  size_t size_of_subsample;
  
  // Not stored here. Stored in "main'.
  DataSubsampler* subsampler;

};

#endif
