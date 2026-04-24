#ifndef STOCHASTICSCANMCMC_H
#define STOCHASTICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
  /**
   * @file stochastic_scan_mcmc.h
   * @brief Defines the StochasticScanMCMC class.
   *
   * A stochastic scan MCMC sampler. Implements the MCMC interface using stochastic scan update steps.
   *
   * @namespace ilike
   * @class StochasticScanMCMC
   * @brief A stochastic scan mcmc derived from MCMC.
   */


class StochasticScanMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Default constructor for StochasticScanMCMC.
   */
  StochasticScanMCMC();
  
  StochasticScanMCMC(const std::vector<MCMC*> &moves,
                     const arma::colvec &unnormalised_probabilities_in);
  
  StochasticScanMCMC(MCMCTermination* termination_in,
                     const std::vector<MCMC*> &moves_in,
                     const arma::colvec &unnormalised_probabilities_in);
  
  /**
   * @brief Destructor for StochasticScanMCMC.
   */
  virtual ~StochasticScanMCMC();
  
  /**
   * @brief Copy constructor for StochasticScanMCMC.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  StochasticScanMCMC(const StochasticScanMCMC &another);
  
  /**
   * @brief Assignment operator for StochasticScanMCMC.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  void operator=(const StochasticScanMCMC &another);
  /**
   * @brief Creates a deep copy of this StochasticScanMCMC object.
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
   * @brief Creates a deep copy and returns it as a stochastic_scan_mcmc pointer.
   *
   * @return The result.
   */
  StochasticScanMCMC* stochastic_scan_mcmc_duplicate() const;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
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
  
  /**
   * @brief Returns the duplicate moves.
   *
   * @return The result.
   */
  std::vector<MCMC*> get_duplicate_moves() const;
  
  // Stored here.
  /** @brief The moves. */
  std::vector<MCMC*> moves;
  
  /** @brief The probabilities. */
  arma::colvec probabilities;
  
  /**
   * @brief Copies the state of another StochasticScanMCMC into this object.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  void make_copy(const StochasticScanMCMC &another);
  
};
}

#endif
