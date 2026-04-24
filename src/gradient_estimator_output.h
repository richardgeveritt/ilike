#ifndef GRADIENTESTIMATOROUTPUT_H
#define GRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "particle.h"

namespace ilike
{
  /**
   * @file gradient_estimator_output.h
   * @brief Defines the GradientEstimatorOutput class.
   *
   * Stores the output of a gradient estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class GradientEstimatorOutput
   * @brief The gradient estimator output class.
   */


class GradientEstimatorOutput
{
  
public:
  
  /**
   * @brief Default constructor for GradientEstimatorOutput.
   */
  GradientEstimatorOutput();
  /**
   * @brief Destructor for GradientEstimatorOutput.
   */
  virtual ~GradientEstimatorOutput();
  
  /**
   * @brief Copy constructor for GradientEstimatorOutput.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  GradientEstimatorOutput(const GradientEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for GradientEstimatorOutput.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void operator=(const GradientEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this GradientEstimatorOutput object.
   *
   * @return The result.
   */
  virtual GradientEstimatorOutput* duplicate() const=0;
  
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Index* index,
                                        const Particle &particle)=0;
  
  //virtual arma::mat get_gradient_of_log_at_previous(const std::string &variable,
  //                                                  const Index* index,
  //                                                  const Particle &particle)=0;
  
  /*
   virtual arma::mat get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters_in)=0;
   */
  
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Index* index,
                                                  const Particle &particle)=0;
  
  //virtual arma::mat subsample_get_gradient_of_log_at_previous(const std::string &variable,
  //                                                            const Index* index,
  //                                                            const Particle &particle)=0;
  
  /**
   * @brief Simulates auxiliary variables.
   */
  virtual void simulate_auxiliary_variables()=0;
  
  /*
   virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters_in)=0;
   */
  
protected:
  
  /**
   * @brief Copies the state of another GradientEstimatorOutput into this object.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void make_copy(const GradientEstimatorOutput &another);
  
};
}

#endif
