#ifndef METROPOLISHASTINGSMCMC_H
#define METROPOLISHASTINGSMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"
#include "proposal_kernel.h"

namespace ilike
{
  /**
   * @file metropolis_hastings_mcmc.h
   * @brief Defines the StandardMCMCOutput class.
   *
   * Stores and manages the output produced by StandardMCMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class StandardMCMCOutput
   * @brief The standard mcmc output class.
   */


class StandardMCMCOutput;

class MetropolisHastingsMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Performs the metropolishastingsmcmc operation.
   */
  MetropolisHastingsMCMC();
  
  // Gaussian random walk.
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         const std::string &variable_name_in,
                         const arma::mat &proposal_covariance_in);
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         const std::vector<std::string> &variable_names_in,
                         const std::vector<arma::mat> &proposal_covariances_in);
  
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         ProposalKernel* proposal_in);
  
  MetropolisHastingsMCMC(MCMCTermination* termination_in,
                         ProposalKernel* proposal_in);
  
  /**
   * @brief Performs the ~metropolishastingsmcmc operation.
   */
  virtual ~MetropolisHastingsMCMC();
  
  /**
   * @brief Performs the metropolishastingsmcmc operation.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  MetropolisHastingsMCMC(const MetropolisHastingsMCMC &another);
  
  /**
   * @brief Assignment operator for StandardMCMCOutput.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  void operator=(const MetropolisHastingsMCMC &another);
  /**
   * @brief Creates a deep copy of this StandardMCMCOutput object.
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
   * @brief Creates a deep copy and returns it as a metropolis_hastings_mcmc pointer.
   *
   * @return The result.
   */
  MetropolisHastingsMCMC* metropolis_hastings_mcmc_duplicate() const;
  
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
  
  /*
   EnsembleMember move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle) const;
   
   EnsembleMember move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle,
   const Parameters &conditioned_on_parameters) const;
   
   EnsembleMember subsample_move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle,
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
  ProposalKernel* proposal;
  
  /**
   * @brief Copies the state of another StandardMCMCOutput into this object.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  void make_copy(const MetropolisHastingsMCMC &another);
  
  // stored here
  /** @brief The index. */
  Index* index;
  
};
}

#endif
