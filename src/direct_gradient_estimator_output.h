#ifndef DIRECTGRADIENTESTIMATOROUTPUT_H
#define DIRECTGRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <boost/unordered_map.hpp>

#include "gradient_estimator_output.h"
#include "particle.h"

class DirectGradientEstimator;

class DirectGradientEstimatorOutput : public GradientEstimatorOutput
{

public:

  DirectGradientEstimatorOutput();
  
  DirectGradientEstimatorOutput(DirectGradientEstimator* estimator_in);

  virtual ~DirectGradientEstimatorOutput();

  DirectGradientEstimatorOutput(const DirectGradientEstimatorOutput &another);

  void operator=(const DirectGradientEstimatorOutput &another);
  GradientEstimatorOutput* duplicate() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle);
  
  /*
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters_in);
  */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle);
  
  /*
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters_in);
  */
  
protected:
  
  // stored in the proposal
  DirectGradientEstimator* estimator;

  void make_copy(const DirectGradientEstimatorOutput &another);
  
  boost::unordered_map< std::string, arma::mat> gradients;
  
  boost::unordered_map< std::string, arma::mat> subsample_gradients;

};

#endif
