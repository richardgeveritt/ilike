#ifndef MCMC_H
#define MCMC_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include "kernel.h"
#include "particle.h"
#include "distributions.h"

namespace ilike
{

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
  
  MCMC();
  MCMC(size_t number_of_iterations_in);
  MCMC(MCMCTermination* termination_in);
  virtual ~MCMC();
  
  MCMC(const MCMC &another);
  
  void operator=(const MCMC &another);
  virtual MCMC* mcmc_duplicate() const=0;
  
  MoveOutput* run(RandomNumberGenerator &rng,
                  const Particle &particle) const;

  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            const Particle &particle) const;

  // One application of the kernel to particle, updating it in place. The
  // result is the same as particle = move(rng, particle), but a rejected
  // proposal leaves particle untouched instead of returning a copy of it
  // (a full copy of the parameters and of every factor output). The default
  // does exactly the assignment; Metropolis, Metropolis-Hastings and the
  // sweeps override it. Used by run().
  virtual void move_in_place(RandomNumberGenerator &rng,
                             Particle &particle) const;

  // As move_in_place(rng, particle), with a scratch particle owned by the
  // caller (one per chain) in which proposals are built and which is
  // swapped with particle on acceptance (dev/plan.md, Stage 3 (a)). Same
  // random numbers and results. Default: ignores scratch.
  virtual void move_in_place(RandomNumberGenerator &rng,
                             Particle &particle,
                             Particle &scratch) const;

  // As move_in_place(rng, particle, scratch), then one step of adaptation of
  // this chain's own copy of the move (dev/plan.md 3.2): used by run() only
  // when the move adapts. Metropolis and Metropolis-Hastings moves give the
  // adaptor the previous, proposed and current states and the acceptance
  // probability; the default refuses (MCMCAdaptor::mcmc_adapt).
  virtual void adaptive_move_in_place(RandomNumberGenerator &rng,
                                      Particle &particle,
                                      Particle &scratch);

  // The adaptation state of this move's kernels (MCMCAdaptor::append_state);
  // nothing for moves that do not adapt.
  virtual void append_adaptation_state(std::vector<std::string> *names,
                                       std::vector<double> *values,
                                       bool final,
                                       const std::string &prefix) const;

  // Per-iteration diagnostics (3.18 item 9, C212): the names of the columns
  // (none by default; nuts: tree_depth, leapfrog_steps, divergent,
  // accept_stat). When there are some, run() moves with
  // diagnostic_move_in_place, which is move_in_place (or, if adapting,
  // adaptive_move_in_place) also filling one row; the row goes to the chain's
  // output, not to the move, which several chains may share.
  virtual void diagnostic_names(std::vector<std::string> &names) const {}
  virtual void diagnostic_move_in_place(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        Particle &scratch,
                                        bool adapting,
                                        double* row);
  // Whether the move is in a warm-up: adapting with `until`, not yet stopped
  // (diagnostics after it count towards the divergence warning).
  virtual bool warming_up() const { return false; }

  // Store every thin-th iteration of a run (MCMC(..., thin); default 1).
  void set_thin(size_t thin_in);
  size_t get_thin() const { return this->thin; }
  // (dev/plan.md 3.21 (c5)) whether a run stores its iterations: an SMC's
  // moves need only each chain's final state (kept apart), so the SMC turns
  // this off (not for MCMC runs, whose chains are their output, nor for
  // waste-free SMC, which uses every state; a chain that adapts or records
  // diagnostics still records)
  void set_store_iterations(bool store_iterations_in) { this->store_iterations = store_iterations_in; }

  // The number of iterations a run makes, when the termination fixes it in
  // advance (iterations(n)); 0 otherwise.
  size_t known_number_of_iterations() const;

  // Whether the move adapts its kernel as it runs (a proposal with an MCMC
  // adaptor). If not, the particles' outputs share one kernel.
  virtual bool adapts() const { return false; }
  
  // could change to also have MoveOutput*
  // change when needed
  // EnsembleMember and particle should merge
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
  virtual void set_index(Index* index_in)=0;
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
  
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);
  
  MCMCTermination* get_duplicated_termination() const;
  MCMCTermination* get_termination() const { return this->termination.get(); }

protected:

  // Every thin-th iteration of a run is stored.
  size_t thin = 1;
  bool store_iterations = true;
  
  friend StochasticScanMCMC;
  friend DeterministicScanMCMC;
  
  virtual StandardMCMCOutput* initialise_mcmc_output() const=0;
  
  virtual void specific_mcmc_adapt(const Particle &current_particle,
                                   size_t iteration_counter)=0;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<MCMCTermination> termination;

  // Owned copy of the cut factors (NULL: none).
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Index> cut_index;

  // Owned union of the move's index and cut_index (NULL when there is no cut):
  // the factors that get new outputs in a proposed particle.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Index> index_with_cut;

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
  
  void make_copy(const MCMC &another);
  
};
}

#endif
