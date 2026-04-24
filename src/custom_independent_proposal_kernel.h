#ifndef CUSTOMINDEPENDENTPROPOSALKERNEL_H
#define CUSTOMINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file custom_independent_proposal_kernel.h
   * @brief Defines the CustomIndependentProposalKernel class.
   *
   * An independent custom proposal kernel. Proposes new parameter values independently of the current state by sampling from a custom distribution.
   *
   * @namespace ilike
   * @class CustomIndependentProposalKernel
   * @brief A custom independent proposal kernel derived from IndependentProposalKernel.
   */


class CustomIndependentProposalKernel : public IndependentProposalKernel
{

public:

  /**
   * @brief Default constructor for CustomIndependentProposalKernel.
   */
  CustomIndependentProposalKernel();
  /**
   * @brief Destructor for CustomIndependentProposalKernel.
   */
  virtual ~CustomIndependentProposalKernel();
  
  /**
   * @brief Constructs a CustomIndependentProposalKernel object.
   *
   * @param proposal_simulate_in The proposal simulate.
   */
  CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in);
  
  CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
                                  EvaluateLogIndependentProposalPtr proposal_evaluate_in);
  
  //CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
  //                                EvaluateLogIndependentProposalPtr proposal_evaluate_in,
  //                                Parameters* proposal_parameters_in);

  /**
   * @brief Copy constructor for CustomIndependentProposalKernel.
   *
   * @param another The CustomIndependentProposalKernel instance to copy from.
   */
  CustomIndependentProposalKernel(const CustomIndependentProposalKernel &another);

  /**
   * @brief Assignment operator for CustomIndependentProposalKernel.
   *
   * @param another The CustomIndependentProposalKernel instance to copy from.
   */
  void operator=(const CustomIndependentProposalKernel &another);
  /**
   * @brief Creates a deep copy of this CustomIndependentProposalKernel object.
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
  //arma::mat independent_gradient_of_log(const std::string &variable,
  //                                      Variables* proposed_particle,
  //                                      const Parameters &conditioned_on_parameters);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  //arma::mat subsample_independent_gradient_of_log(const std::string &variable,
  //                                                Variables* proposed_particle,
  //                                                const Parameters &conditioned_on_parameters);
  
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
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  /**
   * @brief Copies the state of another CustomIndependentProposalKernel into this object.
   *
   * @param another The CustomIndependentProposalKernel instance to copy from.
   */
  void make_copy(const CustomIndependentProposalKernel &another);
  
  /** @brief The proposal evaluate. */
  EvaluateLogIndependentProposalPtr proposal_evaluate;
  
  /** @brief The proposal simulate. */
  SimulateIndependentProposalPtr proposal_simulate;
  
  // not stored here - points to algorithm_parameters
  /** @brief The proposal parameters. */
  Parameters* proposal_parameters;
  
};
}

#endif
