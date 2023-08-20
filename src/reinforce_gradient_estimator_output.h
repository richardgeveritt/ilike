#ifndef REINFORCEGRADIENTESTIMATOROUTPUT_H
#define REINFORCEGRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator_output.h"
#include "particle.h"

class ReinforceGradientEstimator;

class ReinforceInfo
{
public:
  
  ReinforceInfo();
  virtual ~ReinforceInfo();
  
  ReinforceInfo(const ReinforceInfo &another);
  
  void operator=(const ReinforceInfo &another);
  
  // Can populate with all of the auxiliary variables if we need to.
  
  // For now just store the arma::mat gradient;
  arma::mat gradient;
  
protected:
  
  void make_copy(const ReinforceInfo &another);
  
};

class ReinforceGradientEstimatorOutput : public GradientEstimatorOutput
{

public:

  ReinforceGradientEstimatorOutput();
  ReinforceGradientEstimatorOutput(ReinforceGradientEstimator* estimator_in);
  virtual ~ReinforceGradientEstimatorOutput();

  ReinforceGradientEstimatorOutput(const ReinforceGradientEstimatorOutput &another);

  void operator=(const ReinforceGradientEstimatorOutput &another);
  GradientEstimatorOutput* duplicate() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                const Particle &particle);
  
  /*
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters_in);
  */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          const Particle &particle);
  
  void simulate_auxiliary_variables();
  
  /*
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters_in);
  */
  
protected:
  
  // stored in the proposal
  ReinforceGradientEstimator* estimator;
  
  boost::unordered_map< std::string, ReinforceInfo> infos;
  boost::unordered_map< std::string, ReinforceInfo> subsample_infos;
  
  boost::unordered_map< std::string, std::vector<arma::mat>> auxiliary_variables;

  void make_copy(const ReinforceGradientEstimatorOutput &another);

};

#endif
