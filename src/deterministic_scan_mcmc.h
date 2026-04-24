#ifndef DETERMINISTICSCANMCMC_H
#define DETERMINISTICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
  /**
   * @file deterministic_scan_mcmc.h
   * @brief Defines the DeterministicScanMCMC class.
   *
   * A deterministic scan MCMC sampler. Implements the MCMC interface using deterministic scan update steps.
   *
   * @namespace ilike
   * @class DeterministicScanMCMC
   * @brief A deterministic scan mcmc derived from MCMC.
   */


class DeterministicScanMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Default constructor for DeterministicScanMCMC.
   */
  DeterministicScanMCMC();
  
  /**
   * @brief Constructs a DeterministicScanMCMC object.
   *
   * @param moves_in The moves.
   */
  DeterministicScanMCMC(const std::vector<MCMC*> &moves_in);
  
  DeterministicScanMCMC(MCMCTermination* termination_in,
                        const std::vector<MCMC*> &moves_in);
  
  /**
   * @brief Destructor for DeterministicScanMCMC.
   */
  virtual ~DeterministicScanMCMC();
  
  /**
   * @brief Copy constructor for DeterministicScanMCMC.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  DeterministicScanMCMC(const DeterministicScanMCMC &another);
  
  /**
   * @brief Assignment operator for DeterministicScanMCMC.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  void operator=(const DeterministicScanMCMC &another);
  /**
   * @brief Creates a deep copy of this DeterministicScanMCMC object.
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
   * @brief Creates a deep copy and returns it as a deterministic_scan_mcmc pointer.
   *
   * @return The result.
   */
  DeterministicScanMCMC* deterministic_scan_mcmc_duplicate() const;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   
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
  
  /**
   * @brief Returns the duplicate moves.
   *
   * @return The result.
   */
  std::vector<MCMC*> get_duplicate_moves() const;
  
  // Stored here.
  /** @brief The moves. */
  std::vector<MCMC*> moves;
  
  /** @brief The order. */
  std::vector<size_t> order;
  
  /**
   * @brief Copies the state of another DeterministicScanMCMC into this object.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  void make_copy(const DeterministicScanMCMC &another);
  
};
}

#endif
