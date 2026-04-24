#ifndef REINFORCEGRADIENTESTIMATOROUTPUT_H
#define REINFORCEGRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator_output.h"
#include "particle.h"

namespace ilike
{
  /**
   * @file reinforce_gradient_estimator_output.h
   * @brief Defines the ReinforceGradientEstimator class.
   *
   * Estimates the reinforce gradient for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class ReinforceGradientEstimator
   * @brief The reinforce gradient estimator class.
   */


class ReinforceGradientEstimator;

class ReinforceInfo
{
public:
  
  /**
   * @brief Performs the reinforceinfo operation.
   */
  ReinforceInfo();
  /**
   * @brief Performs the ~reinforceinfo operation.
   */
  virtual ~ReinforceInfo();
  
  /**
   * @brief Performs the reinforceinfo operation.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  ReinforceInfo(const ReinforceInfo &another);
  
  /**
   * @brief Assignment operator for ReinforceGradientEstimator.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void operator=(const ReinforceInfo &another);
  
  // Can populate with all of the auxiliary variables if we need to.
  
  // For now just store the arma::mat gradient;
  arma::mat gradient;
  
protected:
  
  /**
   * @brief Copies the state of another ReinforceGradientEstimator into this object.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void make_copy(const ReinforceInfo &another);
  
};

class ReinforceGradientEstimatorOutput : public GradientEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the reinforcegradientestimatoroutput operation.
   */
  ReinforceGradientEstimatorOutput();
  /**
   * @brief Performs the reinforcegradientestimatoroutput operation.
   *
   * @param estimator_in The estimator.
   */
  ReinforceGradientEstimatorOutput(ReinforceGradientEstimator* estimator_in);
  /**
   * @brief Performs the ~reinforcegradientestimatoroutput operation.
   */
  virtual ~ReinforceGradientEstimatorOutput();
  
  /**
   * @brief Performs the reinforcegradientestimatoroutput operation.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  ReinforceGradientEstimatorOutput(const ReinforceGradientEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for ReinforceGradientEstimator.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void operator=(const ReinforceGradientEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this ReinforceGradientEstimator object.
   *
   * @return The result.
   */
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
  
  /**
   * @brief Simulates auxiliary variables.
   */
  void simulate_auxiliary_variables();
  
  /*
   arma::mat subsample_get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters_in);
   */
  
protected:
  
  // stored in the proposal
  /** @brief The estimator. */
  ReinforceGradientEstimator* estimator;
  
  boost::unordered_map< std::string, ReinforceInfo> infos;
  boost::unordered_map< std::string, ReinforceInfo> subsample_infos;
  
  boost::unordered_map< std::string, std::vector<arma::mat>> auxiliary_variables;
  
  /**
   * @brief Copies the state of another ReinforceGradientEstimator into this object.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void make_copy(const ReinforceGradientEstimatorOutput &another);
  
};
}

#endif
