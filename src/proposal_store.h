#ifndef PROPOSALSTORE_H
#define PROPOSALSTORE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "parameters.h"

namespace ilike
{
  /**
   * @file proposal_store.h
   * @brief Defines the GradientEstimatorOutput class.
   *
   * Stores the output of a gradient estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class GradientEstimatorOutput
   * @brief The gradient estimator output class.
   */


class GradientEstimatorOutput;

class ProposalStore
{
  
public:
  
  /**
   * @brief Performs the proposalstore operation.
   */
  ProposalStore();
  /**
   * @brief Performs the ~proposalstore operation.
   */
  virtual ~ProposalStore();
  
  /**
   * @brief Performs the proposalstore operation.
   *
   * @param transformed_parameters_in The transformed parameters.
   */
  ProposalStore(const Parameters &transformed_parameters_in);
  
  /**
   * @brief Performs the proposalstore operation.
   *
   * @param gradient_estimator_output_in The gradient estimator output.
   */
  ProposalStore(GradientEstimatorOutput* gradient_estimator_output_in);
  
  ProposalStore(const Parameters &transformed_parameters_in,
                GradientEstimatorOutput* gradient_estimator_output_in);
  
  /**
   * @brief Performs the proposalstore operation.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  ProposalStore(const ProposalStore &another);
  
  /**
   * @brief Assignment operator for GradientEstimatorOutput.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void operator=(const ProposalStore &another);
  
  /**
   * @brief Sets the transformed parameters.
   *
   * @param transformed_parameters_in The transformed parameters.
   */
  void set_transformed_parameters(const Parameters &transformed_parameters_in);
  /**
   * @brief Sets the gradient estimator output.
   *
   * @param gradient_estimator_output_in The gradient estimator output.
   */
  void set_gradient_estimator_output(GradientEstimatorOutput* gradient_estimator_output_in);
  
  /**
   * @brief Returns the gradient estimator output.
   *
   * @return The result.
   */
  GradientEstimatorOutput* get_gradient_estimator_output() const;
  
  /**
   * @brief Returns the transformed parameters.
   *
   * @return The result.
   */
  Parameters get_transformed_parameters() const;
  
protected:
  
  /** @brief The transformed parameters. */
  Parameters transformed_parameters;
  
  // stored here
  /** @brief The gradient estimator output. */
  GradientEstimatorOutput* gradient_estimator_output;
  
  /**
   * @brief Copies the state of another GradientEstimatorOutput into this object.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void make_copy(const ProposalStore &another);
  
};
}

#endif
