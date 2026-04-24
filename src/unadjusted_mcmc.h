#ifndef UNADJUSTEDMCMC_H
#define UNADJUSTEDMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
  /**
   * @file unadjusted_mcmc.h
   * @brief Defines the EnsembleKalmanOutput class.
   *
   * Stores and manages the output produced by EnsembleKalman. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class EnsembleKalmanOutput
   * @brief The ensemble kalman output class.
   */


class EnsembleKalmanOutput;
class StandardMCMCOutput;
class ProposalKernel;

class UnadjustedMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Performs the unadjustedmcmc operation.
   */
  UnadjustedMCMC();
  
  // Gaussian random walk.
  UnadjustedMCMC(size_t number_of_iterations_in,
                 const std::vector<Parameters> &initial_points_in,
                 const Parameters &proposal_variances);
  
  UnadjustedMCMC(size_t number_of_iterations_in,
                 ProposalKernel* proposal_in);
  
  UnadjustedMCMC(MCMCTermination* termination_in,
                 ProposalKernel* proposal_in);
  
  /**
   * @brief Performs the ~unadjustedmcmc operation.
   */
  virtual ~UnadjustedMCMC();
  
  /**
   * @brief Performs the unadjustedmcmc operation.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  UnadjustedMCMC(const UnadjustedMCMC &another);
  
  /**
   * @brief Assignment operator for EnsembleKalmanOutput.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void operator=(const UnadjustedMCMC &another);
  /**
   * @brief Creates a deep copy of this EnsembleKalmanOutput object.
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
   * @brief Creates a deep copy and returns it as a unadjusted_mcmc pointer.
   *
   * @return The result.
   */
  UnadjustedMCMC* unadjusted_mcmc_duplicate() const;
  
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
  ProposalKernel* proposal;
  
  /**
   * @brief Copies the state of another EnsembleKalmanOutput into this object.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void make_copy(const UnadjustedMCMC &another);
  
  // stored here (change to shared pointer so that memory is not duplicated)
  /** @brief The index. */
  Index* index;
  
};
}

#endif
