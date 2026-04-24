#ifndef ZIGZAGMCMC_H
#define ZIGZAGMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"
#include "proposal_kernel.h"

namespace ilike
{
  /**
   * @file zig_zag_mcmc.h
   * @brief Defines the StandardMCMCOutput class.
   *
   * Stores and manages the output produced by StandardMCMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class StandardMCMCOutput
   * @brief The standard mcmc output class.
   */


// Will not be needed later...
class StandardMCMCOutput;

class ZigZagMCMC : public MCMC
{
  
public:
  
  /**
   * @brief Performs the zigzagmcmc operation.
   */
  ZigZagMCMC();
  
  // Gaussian random walk.
  ZigZagMCMC(size_t number_of_iterations_in,
             const std::vector<Parameters> &initial_points_in,
             const Parameters &proposal_variances);
  
  ZigZagMCMC(size_t number_of_iterations_in,
             ProposalKernel* proposal_in);
  
  /**
   * @brief Performs the ~zigzagmcmc operation.
   */
  virtual ~ZigZagMCMC();
  
  /**
   * @brief Performs the zigzagmcmc operation.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  ZigZagMCMC(const ZigZagMCMC &another);
  
  /**
   * @brief Assignment operator for StandardMCMCOutput.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  void operator=(const ZigZagMCMC &another);
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
  
  // Zig zag can be run to produce a Particle, just as any MCMC can be. This way is can be combined with other MCMC moves.
  // When we call "run" on it, as in the base class, it will use a succession of moves and produce a succession of particles, just as any other MCMC move would.
  // Instead of this, we need to introduce a base class for standard MCMC and another for PDMPs.
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
  
  // Overrides that just run the PDMP, rather than calling run in the MCMC base class.
  MoveOutput* run(RandomNumberGenerator &rng,
                  const Particle &particle) const; // override?
  
  /*
   MoveOutput* run(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const; // override?
   */
  
  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            const Particle &particle) const;
  
  /*
   MoveOutput* subsample_run(RandomNumberGenerator &rng,
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
   * @brief Copies the state of another StandardMCMCOutput into this object.
   *
   * @param another The StandardMCMCOutput instance to copy from.
   */
  void make_copy(const ZigZagMCMC &another);
  
};
}

#endif
