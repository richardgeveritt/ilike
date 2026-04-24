#ifndef DIRECTGRADIENTESTIMATOROUTPUT_H
#define DIRECTGRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <boost/unordered_map.hpp>

#include "gradient_estimator_output.h"
#include "particle.h"

namespace ilike
{
  /**
   * @file direct_gradient_estimator_output.h
   * @brief Defines the DirectGradientEstimator class.
   *
   * Estimates the direct gradient for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class DirectGradientEstimator
   * @brief The direct gradient estimator class.
   */


class DirectGradientEstimator;

class DirectGradientEstimatorOutput : public GradientEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the directgradientestimatoroutput operation.
   */
  DirectGradientEstimatorOutput();
  
  /**
   * @brief Performs the directgradientestimatoroutput operation.
   *
   * @param estimator_in The estimator.
   */
  DirectGradientEstimatorOutput(DirectGradientEstimator* estimator_in);
  
  /**
   * @brief Performs the ~directgradientestimatoroutput operation.
   */
  virtual ~DirectGradientEstimatorOutput();
  
  /**
   * @brief Performs the directgradientestimatoroutput operation.
   *
   * @param another The DirectGradientEstimator instance to copy from.
   */
  DirectGradientEstimatorOutput(const DirectGradientEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for DirectGradientEstimator.
   *
   * @param another The DirectGradientEstimator instance to copy from.
   */
  void operator=(const DirectGradientEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this DirectGradientEstimator object.
   *
   * @return The result.
   */
  GradientEstimatorOutput* duplicate() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                const Particle &particle);
  
  //arma::mat get_gradient_of_log_at_previous(const std::string &variable,
  //                                          const Index* index,
  //                                          const Particle &particle);
  
  /*
   arma::mat get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters_in);
   */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          const Particle &particle);
  
  //arma::mat subsample_get_gradient_of_log_at_previous(const std::string &variable,
  //                                                    const Index* index,
  //                                                    const Particle &particle);
  
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
  DirectGradientEstimator* estimator;
  
  /**
   * @brief Copies the state of another DirectGradientEstimator into this object.
   *
   * @param another The DirectGradientEstimator instance to copy from.
   */
  void make_copy(const DirectGradientEstimatorOutput &another);
  
  boost::unordered_map< std::string, arma::mat> gradients;
  
  boost::unordered_map< std::string, arma::mat> subsample_gradients;
  
};
}

#endif
