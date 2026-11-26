#ifndef RATIOESTIMATINGMHMCMC_H
#define RATIOESTIMATINGMHMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "metropolis_hastings_mcmc.h"
#include "estimated_ratio.h"

namespace ilike
{

// As RatioEstimatingMHMCMC, with a proposal that need not be symmetric (its
// densities enter the ratio as in MetropolisHastingsMCMC): for example the
// exchange algorithm with a Langevin proposal whose gradient includes an
// estimate of the gradient of log 1/Z (MALA-exchange, Alquier, Friel, Everitt
// & Boland 2016; dev/plan.md 3.22 U1b, C287). The estimated term is drawn
// independently of the fields behind the gradient: reusing one of those
// fields, as Algorithm 9 of that paper does, does not leave the target
// invariant (checked by enumeration, C287).
class RatioEstimatingMHMCMC : public MetropolisHastingsMCMC
{

public:

  RatioEstimatingMHMCMC();

  RatioEstimatingMHMCMC(size_t number_of_iterations_in,
                      ProposalKernel* proposal_in,
                      RatioStep* step_in);

  RatioEstimatingMHMCMC(MCMCTermination* termination_in,
                      ProposalKernel* proposal_in,
                      RatioStep* step_in);

  virtual ~RatioEstimatingMHMCMC();

  RatioEstimatingMHMCMC(const RatioEstimatingMHMCMC &another);

  void operator=(const RatioEstimatingMHMCMC &another);
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

  // (no per-iteration diagnostics: those of hmc describe its own accept step)
  void diagnostic_names(std::vector<std::string> &names) const { names.clear(); }

  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;

  // points = N: the averaged acceptance ratio (above); parallel: the N draws
  // in RcppParallel's threads (the same results either way)
  EstimatedRatio& estimated_ratio() { return this->ratio; }

protected:

  StandardMCMCOutput* initialise_mcmc_output() const;

  // the decision for a proposed particle (after the uniform is drawn): the
  // move's factors first (a proposal they rule out draws no field), then the
  // fields' terms. Sets the acceptance probability, for tuning and adaptation.
  bool decide(RandomNumberGenerator &rng,
              double log_u,
              const Particle &particle,
              Particle &proposed_particle,
              double &acceptance_probability) const;

  void accept(Particle &proposed_particle) const;

  void make_copy(const RatioEstimatingMHMCMC &another);

  // the estimated term (the step, the number of draws, ...)
  EstimatedRatio ratio;

};
}

#endif
