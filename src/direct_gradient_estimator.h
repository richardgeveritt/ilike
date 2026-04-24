#ifndef DIRECTGRADIENTESTIMATOR_H
#define DIRECTGRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator.h"
#include "particle.h"

namespace ilike
{
  /**
   * @file direct_gradient_estimator.h
   * @brief Defines the GradientEstimatorOutput class.
   *
   * Stores the output of a gradient estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class GradientEstimatorOutput
   * @brief The gradient estimator output class.
   */


class GradientEstimatorOutput;

class DirectGradientEstimator : public GradientEstimator
{
  
public:
  
  /**
   * @brief Performs the directgradientestimator operation.
   */
  DirectGradientEstimator();
  
  /**
   * @brief Performs the ~directgradientestimator operation.
   */
  virtual ~DirectGradientEstimator();
  
  /**
   * @brief Performs the directgradientestimator operation.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  DirectGradientEstimator(const DirectGradientEstimator &another);
  
  /**
   * @brief Assignment operator for GradientEstimatorOutput.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void operator=(const DirectGradientEstimator &another);
  /**
   * @brief Creates a deep copy of this GradientEstimatorOutput object.
   *
   * @return The result.
   */
  GradientEstimator* duplicate() const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  GradientEstimatorOutput* initialise();
  
  //GradientEstimatorOutput* generate_new_gradient_estimator_output(ProposalKernel* proposal,
  //                                                                Particle &particle);
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Index* index,
                                const Particle &particle) const;
  
  /*
   arma::mat get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Index* index,
                                          const Particle &particle) const;
  
  /*
   arma::mat subsample_get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
protected:
  
  /**
   * @brief Copies the state of another GradientEstimatorOutput into this object.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void make_copy(const DirectGradientEstimator &another);
  
};
}

#endif
