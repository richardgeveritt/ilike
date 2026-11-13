#ifndef MCMC_H
#define MCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "kernel.h"
#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file mcmc.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class MCMCTermination;
class MoveOutput;
class MCMCAdaptor;
class StochasticScanMCMC;
class DeterministicScanMCMC;
class Index;
class StandardMCMCOutput;

class MCMC : public Kernel
{
  
public:
  
  /**
   * @brief Performs the mcmc operation.
   */
  MCMC();
  /**
   * @brief Performs the mcmc operation.
   *
   * @param number_of_iterations_in The number of iterations.
   */
  MCMC(size_t number_of_iterations_in);
  /**
   * @brief Performs the mcmc operation.
   *
   * @param termination_in The termination.
   */
  MCMC(MCMCTermination* termination_in);
  /**
   * @brief Performs the ~mcmc operation.
   */
  virtual ~MCMC();
  
  /**
   * @brief Performs the mcmc operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  MCMC(const MCMC &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const MCMC &another);
  /**
   * @brief Creates a deep copy and returns it as a mcmc pointer.
   *
   * @return The result.
   */
  virtual MCMC* mcmc_duplicate() const=0;
  
  MoveOutput* run(RandomNumberGenerator &rng,
                  const Particle &particle) const;
  
  /*
   MoveOutput* run(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            const Particle &particle) const;
  
  /*
   MoveOutput* subsample_run(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  // could change to also have MoveOutput*
  // change when needed
  // EnsembleMember and particle should merge
  /*
   EnsembleMember run(RandomNumberGenerator &rng,
   EnsembleMember &particle);
   
   EnsembleMember run(RandomNumberGenerator &rng,
   EnsembleMember &particle,
   const Parameters &conditioned_on_parameters);
   
   EnsembleMember subsample_run(RandomNumberGenerator &rng,
   EnsembleMember &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
  /**
   * @brief Sets the index.
   *
   * @param index_in The index.
   */
  virtual void set_index(Index* index_in)=0;
  /**
   * @brief Sets the index if null.
   *
   * @param index_in The index.
   */
  virtual void set_index_if_null(Index* index_in)=0;

  // Cut factors: they depend on what this move changes, but are deliberately
  // left out of its acceptance ratio (as in a "cut" model). They are
  // re-evaluated at the new point when a proposal is accepted, so the
  // per-factor cache stays valid for later moves. Copies cut_index_in (the
  // caller keeps ownership). Must be called after the move's own index is set.
  // Only Metropolis, Metropolis-Hastings and unadjusted moves support cuts.
  virtual void set_cut_index(const Index* cut_index_in);

  // Turn on the optional consistency check of the factor indices (see
  // check_carried_over_factors). label names the move in error messages.
  void set_check_factor_index(bool check, const std::string &label);

  // Whether this move reads cached factor values (for an acceptance-ratio
  // denominator). True for Metropolis and Metropolis-Hastings moves; false
  // for unadjusted moves; a sweep reads them if any of its moves does.
  virtual bool reads_cached_factors() const;
  
  /**
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);
  
  /**
   * @brief Returns the duplicated termination.
   *
   * @return The result.
   */
  MCMCTermination* get_duplicated_termination() const;
  
  //virtual void set_termination_parameters(h)=0;
  
protected:
  
  friend StochasticScanMCMC;
  friend DeterministicScanMCMC;
  
  /**
   * @brief Performs the initialise mcmc output operation.
   *
   * @return The result.
   */
  virtual StandardMCMCOutput* initialise_mcmc_output() const=0;
  
  virtual void specific_mcmc_adapt(const Particle &current_particle,
                                   size_t iteration_counter)=0;
  
  // Stored here.
  /** @brief The termination. */
  MCMCTermination* termination;

  // Owned copy of the cut factors (NULL: none).
  Index* cut_index;

  // Owned union of the move's index and cut_index (NULL when there is no cut):
  // the factors that get new outputs in a proposed particle.
  Index* index_with_cut;

  // Recompute index_with_cut from the move's index; call whenever either
  // index changes.
  void update_index_with_cut(const Index* move_index);

  // Replace cut_index by a copy of cut_index_in (used by the subclasses'
  // set_cut_index, which then call update_index_with_cut).
  void store_cut_index(const Index* cut_index_in);

  // The factors that get new outputs in a proposed particle for this move.
  const Index* get_new_outputs_index(const Index* move_index) const;

  // Evaluate the cut factors of an accepted particle at its new parameters.
  void refresh_cut_factors(Particle &accepted_particle) const;

  // Optional check (off by default): after a move is accepted, stop if any
  // exact factor that was carried over (neither in the move's index nor cut)
  // changed value, i.e. it did not cancel in the acceptance ratio.
  bool check_factor_index;
  std::string move_label;
  void check_carried_over_factors(const Particle &accepted_particle,
                                  const Index* move_index) const;
  
  // stored here
  //Index* index;
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const MCMC &another);
  
};
}

#endif
