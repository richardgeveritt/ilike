#ifndef CUSTOMGUIDEDNOPARAMSSYMMETRICPROPOSALKERNEL_H
#define CUSTOMGUIDEDNOPARAMSSYMMETRICPROPOSALKERNEL_H

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
   * @file custom_guided_no_params_symmetric_proposal_kernel.h
   * @brief Defines the CustomGuidedNoParamsSymmetricProposalKernel class.
   *
   * A custom guided no params symmetric proposal kernel. Proposes new parameter values during MCMC or SMC moves using a custom guided no params symmetric distribution centred on the current state.
   *
   * @namespace ilike
   * @class CustomGuidedNoParamsSymmetricProposalKernel
   * @brief A custom guided no params symmetric proposal kernel derived from SymmetricProposalKernel.
   */


class CustomGuidedNoParamsSymmetricProposalKernel : public SymmetricProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for CustomGuidedNoParamsSymmetricProposalKernel.
   */
  CustomGuidedNoParamsSymmetricProposalKernel();
  /**
   * @brief Destructor for CustomGuidedNoParamsSymmetricProposalKernel.
   */
  virtual ~CustomGuidedNoParamsSymmetricProposalKernel();
  
  CustomGuidedNoParamsSymmetricProposalKernel(SimulateGuidedNoParamsMCMCProposalPtr proposal_simulate_in,
                                              Data* data_in);
  
  /**
   * @brief Copy constructor for CustomGuidedNoParamsSymmetricProposalKernel.
   *
   * @param another The CustomGuidedNoParamsSymmetricProposalKernel instance to copy from.
   */
  CustomGuidedNoParamsSymmetricProposalKernel(const CustomGuidedNoParamsSymmetricProposalKernel &another);
  
  /**
   * @brief Assignment operator for CustomGuidedNoParamsSymmetricProposalKernel.
   *
   * @param another The CustomGuidedNoParamsSymmetricProposalKernel instance to copy from.
   */
  void operator=(const CustomGuidedNoParamsSymmetricProposalKernel &another);
  /**
   * @brief Creates a deep copy of this CustomGuidedNoParamsSymmetricProposalKernel object.
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
  
  /*
   double specific_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  
  /*
   double specific_subsample_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  
  /*
   Parameters simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  
  /*
   Parameters subsample_simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  
  /*
   Parameters subsample_simulate(RandomNumberGenerator &rng,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  /*
   arma::mat specific_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;
  
  /*
   arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void make_copy(const CustomGuidedNoParamsSymmetricProposalKernel &another);
  
  /** @brief The proposal evaluate. */
  EvaluateLogGuidedNoParamsMCMCProposalPtr proposal_evaluate;
  
  /** @brief The proposal simulate. */
  SimulateGuidedNoParamsMCMCProposalPtr proposal_simulate;
  
  /** @brief The data. */
  Data* data;
  
};
}

#endif
