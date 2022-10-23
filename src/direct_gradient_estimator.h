#ifndef DIRECTGRADIENTESTIMATOR_H
#define DIRECTGRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator.h"
#include "particle.h"

class GradientEstimatorOutput;

class DirectGradientEstimator : public GradientEstimator
{

public:

  DirectGradientEstimator();

  virtual ~DirectGradientEstimator();

  DirectGradientEstimator(const DirectGradientEstimator &another);

  void operator=(const DirectGradientEstimator &another);
  GradientEstimator* duplicate() const;
  
  GradientEstimatorOutput* initialise() const;
  
  //GradientEstimatorOutput* generate_new_gradient_estimator_output(ProposalKernel* proposal,
  //                                                                Particle &particle);
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle);
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters);
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle);
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters);
  
protected:

  void make_copy(const DirectGradientEstimator &another);

};

#endif
