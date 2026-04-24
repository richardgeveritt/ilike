#ifndef METROPOLISMCMC_H
#define METROPOLISMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
  /**
   * @file metropolis_mcmc.h
   * @brief Defines the SymmetricProposalKernel class.
   *
   * A symmetric proposal kernel. Proposes new parameter values during MCMC or SMC moves using a symmetric distribution centred on the current state.
   *
   * @namespace ilike
   * @class SymmetricProposalKernel
   * @brief The symmetric proposal kernel class.
   */


class SymmetricProposalKernel;
class EnsembleKalmanOutput;
class StandardMCMCOutput;

class MetropolisMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Performs the metropolismcmc operation.
   */
  MetropolisMCMC();
  
  // Gaussian random walk.
  MetropolisMCMC(size_t number_of_iterations_in,
                 const std::vector<Parameters> &initial_points_in,
                 const Parameters &proposal_variances);
  
  MetropolisMCMC(size_t number_of_iterations_in,
                 SymmetricProposalKernel* proposal_in);
  
  MetropolisMCMC(MCMCTermination* termination_in,
                 SymmetricProposalKernel* proposal_in);
  
  /**
   * @brief Performs the ~metropolismcmc operation.
   */
  virtual ~MetropolisMCMC();
  
  /**
   * @brief Performs the metropolismcmc operation.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  MetropolisMCMC(const MetropolisMCMC &another);
  
  /**
   * @brief Assignment operator for SymmetricProposalKernel.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  void operator=(const MetropolisMCMC &another);
  /**
   * @brief Creates a deep copy of this SymmetricProposalKernel object.
   *
   * @return The result.
   */
  Kernel* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a mcmc pointer.
   *
   * @return The result.
   */
  MCMC* mcmc_duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a metropolis_mcmc pointer.
   *
   * @return The result.
   */
  MetropolisMCMC* metropolis_mcmc_duplicate() const;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  void smc_adapt(SMCOutput* current_state);
  /**
   * @brief Performs the ensemble adapt operation.
   *
   * @param current_state The current state.
   */
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
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
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  std::vector<const ProposalKernel*> get_proposals() const;
  
protected:
  
  void specific_mcmc_adapt(const Particle &current_particle,
                           size_t iteration_counter);
  
  /**
   * @brief Performs the initialise mcmc output operation.
   *
   * @return The result.
   */
  StandardMCMCOutput* initialise_mcmc_output() const;
  
  // stored here
  /** @brief The proposal. */
  SymmetricProposalKernel* proposal;
  
  /**
   * @brief Copies the state of another SymmetricProposalKernel into this object.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  void make_copy(const MetropolisMCMC &another);
  
  // stored here (change to shared pointer so that memory is not duplicated)
  /** @brief The index. */
  Index* index;
  
};
}

#endif
