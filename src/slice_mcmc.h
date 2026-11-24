#ifndef SLICEMCMC_H
#define SLICEMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "metropolis_hastings_mcmc.h"

namespace ilike
{

class SlicePointProposalKernel;

// move m = slice(x, width = w, lower = a, upper = b, max_steps = m)
// (dev/plan.md 3.21 (c3); Neal 2003, univariate slice sampling): each element
// of x in turn is moved by a slice sampling update of the move's target (its
// factors, as for any move): a level log y = log pi(x) + log U under the
// current point, an interval around the element — the whole of [a, b] when
// no width is given (both limits finite: e.g. u on [0, 1], as in Prangle,
// Everitt & Kypraios's rare-event ABC), otherwise stepped out by w (at most m
// steps, Neal's procedure) within [a, b] — and points drawn uniformly in it,
// shrinking it towards the current element, until one is on the slice. Each
// update leaves the target invariant, with no rejection. The target must be
// exact (not an estimate). Lowered as a Metropolis-Hastings move whose
// proposal holds the settings (SlicePointProposalKernel), as nuts is.
class SliceMCMC : public MetropolisHastingsMCMC
{
public:

  SliceMCMC();
  SliceMCMC(size_t number_of_iterations_in, SlicePointProposalKernel* proposal_in);
  SliceMCMC(MCMCTermination* termination_in, SlicePointProposalKernel* proposal_in);
  virtual ~SliceMCMC();

  SliceMCMC(const SliceMCMC &another);
  void operator=(const SliceMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  MetropolisHastingsMCMC* metropolis_hastings_mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng,
                     Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng,
                     Particle &particle,
                     Particle &scratch) const;
  void adaptive_move_in_place(RandomNumberGenerator &rng,
                              Particle &particle,
                              Particle &scratch);
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
};

}

#endif
