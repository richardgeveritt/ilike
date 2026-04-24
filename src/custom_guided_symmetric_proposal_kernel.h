#ifndef CUSTOMGUIDEDSYMMETRICPROPOSALKERNEL_H
#define CUSTOMGUIDEDSYMMETRICPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file custom_guided_symmetric_proposal_kernel.h
   * @brief Defines the CustomGuidedSymmetricProposalKernel class.
   *
   * A custom guided symmetric proposal kernel. Proposes new parameter values during MCMC or SMC moves using a custom guided symmetric distribution centred on the current state.
   *
   * @namespace ilike
   * @class CustomGuidedSymmetricProposalKernel
   * @brief A custom guided symmetric proposal kernel derived from SymmetricProposalKernel.
   */


class CustomGuidedSymmetricProposalKernel : public SymmetricProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for CustomGuidedSymmetricProposalKernel.
   */
  CustomGuidedSymmetricProposalKernel();
  /**
   * @brief Destructor for CustomGuidedSymmetricProposalKernel.
   */
  virtual ~CustomGuidedSymmetricProposalKernel();
  
  CustomGuidedSymmetricProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                                      Data* data_in);
  
  /**
   * @brief Copy constructor for CustomGuidedSymmetricProposalKernel.
   *
   * @param another The CustomGuidedSymmetricProposalKernel instance to copy from.
   */
  CustomGuidedSymmetricProposalKernel(const CustomGuidedSymmetricProposalKernel &another);
  
  /**
   * @brief Assignment operator for CustomGuidedSymmetricProposalKernel.
   *
   * @param another The CustomGuidedSymmetricProposalKernel instance to copy from.
   */
  void operator=(const CustomGuidedSymmetricProposalKernel &another);
  /**
   * @brief Creates a deep copy of this CustomGuidedSymmetricProposalKernel object.
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
   * @brief Creates a deep copy and returns it as a symmetric_proposal_kernel pointer.
   *
   * @return The result.
   */
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
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
  
  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  
  
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;
  
  /**
   * @brief Copies the state of another CustomGuidedSymmetricProposalKernel into this object.
   *
   * @param another The CustomGuidedSymmetricProposalKernel instance to copy from.
   */
  void make_copy(const CustomGuidedSymmetricProposalKernel &another);
  
  /** @brief The proposal evaluate. */
  EvaluateLogGuidedMCMCProposalPtr proposal_evaluate;
  
  /** @brief The proposal simulate. */
  SimulateGuidedMCMCProposalPtr proposal_simulate;
  
  /** @brief The proposal parameters. */
  Parameters* proposal_parameters;
  
  /** @brief The data. */
  Data* data;
};
}

#endif
