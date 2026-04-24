#ifndef CUSTOMNOPARAMSPROPOSALKERNEL_H
#define CUSTOMNOPARAMSPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file custom_no_params_proposal_kernel.h
   * @brief Defines the CustomNoParamsProposalKernel class.
   *
   * A custom no params proposal kernel. Proposes new parameter values during MCMC or SMC moves using a custom no params distribution centred on the current state.
   *
   * @namespace ilike
   * @class CustomNoParamsProposalKernel
   * @brief A custom no params proposal kernel derived from ProposalKernel.
   */


class CustomNoParamsProposalKernel : public ProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for CustomNoParamsProposalKernel.
   */
  CustomNoParamsProposalKernel();
  /**
   * @brief Destructor for CustomNoParamsProposalKernel.
   */
  virtual ~CustomNoParamsProposalKernel();
  
  /**
   * @brief Constructs a CustomNoParamsProposalKernel object.
   *
   * @param proposal_simulate_in The proposal simulate.
   */
  CustomNoParamsProposalKernel(SimulateNoParamsMCMCProposalPtr proposal_simulate_in);
  
  CustomNoParamsProposalKernel(SimulateNoParamsMCMCProposalPtr proposal_simulate_in,
                               EvaluateLogNoParamsMCMCProposalPtr proposal_evaluate_in);
  
  /**
   * @brief Copy constructor for CustomNoParamsProposalKernel.
   *
   * @param another The CustomNoParamsProposalKernel instance to copy from.
   */
  CustomNoParamsProposalKernel(const CustomNoParamsProposalKernel &another);
  
  /**
   * @brief Assignment operator for CustomNoParamsProposalKernel.
   *
   * @param another The CustomNoParamsProposalKernel instance to copy from.
   */
  void operator=(const CustomNoParamsProposalKernel &another);
  /**
   * @brief Creates a deep copy of this CustomNoParamsProposalKernel object.
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
  
  void make_copy(const CustomNoParamsProposalKernel &another);
  
  /** @brief The proposal evaluate. */
  EvaluateLogNoParamsMCMCProposalPtr proposal_evaluate;
  
  /** @brief The proposal simulate. */
  SimulateNoParamsMCMCProposalPtr proposal_simulate;
  
};
}

#endif
