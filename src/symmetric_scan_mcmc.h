#ifndef SYMMETRICSCANMCMC_H
#define SYMMETRICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "deterministic_scan_mcmc.h"

namespace ilike
{

// A symmetric scan of moves (dev/plan.md 3.22, "Every annealed path on the
// general engine", stage A): each iteration runs the moves in order and then
// in reverse order, each reversed (MCMC::move_in_place_reversed), which makes
// one iteration a reversible kernel when each move has reversible parts. A
// path inside a Metropolis-Hastings move (AIS MCMC, annealed reversible jump,
// annealed exchange) is exact only with reversible kernels and kernels that
// do not change during the path, so this scan does not adapt: its moves keep
// their settings (adapts() is false, and adaptation calls move without
// adapting). Made by make_mcmc for a recipe with symmetric_scan = true, which
// refuses moves without reversible parts.
//
// (A2) Linked to a cycle, it runs the cycle's own moves (not owned) instead of
// its own: the moves a chain adapts, with the values reached so far, fixed
// for the path. Calls that change the moves' state (indices, storing,
// adaptation) still go to its own moves, so the path's SMC cannot disturb the
// cycle. During each move it gives the particle the caller's tuning (a
// particle's own scales under esjd in SMC, which the moves read by their
// instance), restoring the particle's own afterwards.
struct ParticleTuning;

class SymmetricScanMCMC : public DeterministicScanMCMC
{
public:

  SymmetricScanMCMC();
  SymmetricScanMCMC(MCMCTermination* termination_in, const std::vector<MCMC*> &moves_in);
  virtual ~SymmetricScanMCMC();
  SymmetricScanMCMC(const SymmetricScanMCMC &another);

  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  DeterministicScanMCMC* deterministic_scan_mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng, const Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng, Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch) const;
  void adaptive_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch);
  void diagnostic_names(std::vector<std::string> &names) const { names.clear(); }
  void diagnostic_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch,
                                bool adapting, double* row);
  bool adapts() const { return false; }
  void smc_adapt(SMCOutput* current_state) {}
  // (itself symmetric: its reverse is itself)
  void move_in_place_reversed(RandomNumberGenerator &rng, Particle &particle) const;
  Particle subsample_move(RandomNumberGenerator &rng, const Particle &particle) const;

  // (A2) the moves run (NULL entries are not allowed); empty: its own
  void link(const std::vector<const MCMC*> &moves_in);
  bool linked() const { return !this->linked_moves.empty(); }
  // the proposals of the moves it runs (the particles' stores are keyed by them)
  std::vector<const ProposalKernel*> get_proposals() const;
  // (A2) the tuning the particle carries during each move (NULL: its own)
  void set_caller_tuning(const ParticleTuning* tuning_in) { this->caller_tuning = tuning_in; }

protected:

  void forward(RandomNumberGenerator &rng, Particle &particle, Particle* scratch) const;
  void backward(RandomNumberGenerator &rng, Particle &particle) const;

  // not owned: the cycle's moves
  std::vector<const MCMC*> linked_moves;
  const ParticleTuning* caller_tuning = NULL;

  StandardMCMCOutput* initialise_mcmc_output() const;
};

}

#endif
