#ifndef TRANSFORMEDINDEPENDENTPROPOSALKERNEL_H
#define TRANSFORMEDINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "transform.h"

namespace ilike
{
  /**
   * @file transformed_independent_proposal_kernel.h
   * @brief Defines the TransformedIndependentProposalKernel class.
   *
   * An independent transformed proposal kernel. Proposes new parameter values independently of the current state by sampling from a transformed distribution.
   *
   * @namespace ilike
   * @class TransformedIndependentProposalKernel
   * @brief A transformed independent proposal kernel derived from IndependentProposalKernel.
   */


class TransformedIndependentProposalKernel : public IndependentProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for TransformedIndependentProposalKernel.
   */
  TransformedIndependentProposalKernel();
  /**
   * @brief Destructor for TransformedIndependentProposalKernel.
   */
  virtual ~TransformedIndependentProposalKernel();
  
  TransformedIndependentProposalKernel(IndependentProposalKernel* proposal_in,
                                       std::shared_ptr<Transform> transform_in,
                                       bool distribution_on_transformed_space_in=true);
  
  /**
   * @brief Copy constructor for TransformedIndependentProposalKernel.
   *
   * @param another The TransformedIndependentProposalKernel instance to copy from.
   */
  TransformedIndependentProposalKernel(const TransformedIndependentProposalKernel &another);
  
  /**
   * @brief Assignment operator for TransformedIndependentProposalKernel.
   *
   * @param another The TransformedIndependentProposalKernel instance to copy from.
   */
  void operator=(const TransformedIndependentProposalKernel &another);
  /**
   * @brief Creates a deep copy of this TransformedIndependentProposalKernel object.
   *
   * @return The result.
   */
  Kernel* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a proposal_kernel pointer.
   *
   * @return The result.
   */
  ProposalKernel* proposal_kernel_duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a independent_proposal_kernel pointer.
   *
   * @return The result.
   */
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;
  
  /**
   * @brief Evaluates the independent kernel.
   *
   * @param proposed_particle The proposed particle.
   *
   * @return The result.
   */
  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double evaluate_independent_kernel(Variables* proposed_particle,
  //                                   const Parameters &conditioned_on_parameters) const;
  /**
   * @brief Performs the subsample evaluate independent kernel operation.
   *
   * @param proposed_particle The proposed particle.
   *
   * @return The result.
   */
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  
  /**
   * @brief Performs the independent simulate operation.
   *
   * @param rng The rng.
   *
   * @return The result.
   */
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng,
                                  const Parameters &conditioned_on_parameters) const;
  
  /**
   * @brief Performs the subsample independent simulate operation.
   *
   * @param rng The rng.
   *
   * @return The result.
   */
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable,
                                            const Parameters &conditioned_on_parameters) const;
  
  arma::mat independent_gradient_of_log(const std::string &variable,
                                        const Parameters &proposed_particle);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  
  /**
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  /**
   * @brief Simulates gradient estimator output.
   *
   * @return The result.
   */
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  
  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  std::vector<const ProposalKernel*> get_proposals() const;
  
  /**
   * @brief Sets the index.
   *
   * @param index_in The index.
   */
  void set_index(Index* index_in);
  /**
   * @brief Sets the index if null.
   *
   * @param index_in The index.
   */
  void set_index_if_null(Index* index_in);
  
  /**
   * @brief Performs the can be evaluated operation.
   *
   * @return The result.
   */
  bool can_be_evaluated() const;
  
  /**
   * @brief Sets the data.
   *
   * @param data_in The data.
   */
  void set_data(Data* data_in);
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, t**ake in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  /**
   * @brief Copies the state of another TransformedIndependentProposalKernel into this object.
   *
   * @param another The TransformedIndependentProposalKernel instance to copy from.
   */
  void make_copy(const TransformedIndependentProposalKernel &another);
  
  // stored here
  /** @brief The proposal. */
  IndependentProposalKernel* proposal;
  
  /** @brief The transform. */
  std::shared_ptr<Transform> transform;
  
  /** @brief The distribution on transformed space. */
  bool distribution_on_transformed_space;
  
};
}

#endif
