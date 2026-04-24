#ifndef TRANSFORMEDPROPOSALKERNEL_H
#define TRANSFORMEDPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "transform.h"

namespace ilike
{
  /**
   * @file transformed_proposal_kernel.h
   * @brief Defines the TransformedProposalKernel class.
   *
   * A transformed proposal kernel. Proposes new parameter values during MCMC or SMC moves using a transformed distribution centred on the current state.
   *
   * @namespace ilike
   * @class TransformedProposalKernel
   * @brief A transformed proposal kernel derived from ProposalKernel.
   */


class TransformedProposalKernel : public ProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for TransformedProposalKernel.
   */
  TransformedProposalKernel();
  /**
   * @brief Destructor for TransformedProposalKernel.
   */
  virtual ~TransformedProposalKernel();
  
  TransformedProposalKernel(ProposalKernel* proposal_in,
                            const Transform &transform_in,
                            bool distribution_on_transformed_space_in);
  
  /**
   * @brief Copy constructor for TransformedProposalKernel.
   *
   * @param another The TransformedProposalKernel instance to copy from.
   */
  TransformedProposalKernel(const TransformedProposalKernel &another);
  
  /**
   * @brief Assignment operator for TransformedProposalKernel.
   *
   * @param another The TransformedProposalKernel instance to copy from.
   */
  void operator=(const TransformedProposalKernel &another);
  /**
   * @brief Creates a deep copy of this TransformedProposalKernel object.
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
  
  /*
   arma::mat specific_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  /*
   arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void make_copy(const TransformedProposalKernel &another);
  
  // stored here
  /** @brief The proposal. */
  ProposalKernel* proposal;
  
  /** @brief The transform. */
  Transform transform;
  
  /** @brief The distribution on transformed space. */
  bool distribution_on_transformed_space;
  
};

}

#endif
