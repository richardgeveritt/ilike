#ifndef RATIOESTIMATINGMCMC_H
#define RATIOESTIMATINGMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "metropolis_mcmc.h"
#include "estimated_ratio.h"

namespace ilike
{

// A Metropolis move whose acceptance ratio is estimated (dev/plan.md 3.22 A,
// U1, U2): it proposes vartheta with a symmetric kernel, evaluates the
// move's own factors there (a proposal they rule out draws nothing), and
// multiplies their ratio by the term of an auxiliary draw (a RatioStep; the
// exchange algorithm's is the first). With points = N it uses the averaged
// acceptance ratio of Andrieu, Yildirim, Doucet & Chopin (MHAAR, arXiv
// 2101.01253 Algorithm 1), which keeps the move exact: a direction c, equally
// likely 1 or 2, then
//  c = 1: N draws u_i ~ Q(theta -> vartheta); log r = (factors) + log (1/N) sum_i exp(t_i),
//         t_i = step.draw(theta, vartheta);
//  c = 2: one draw forward, t_1 = -step.draw(theta, vartheta) (the reverse move's
//         term for the same u), and N - 1 backward, t_i = step.draw(vartheta, theta);
//         log r = (factors) - log (1/N) sum_i exp(t_i).
// For N = 1 both directions give the same ratio, so no direction is drawn.
// The N draws of a proposal use their own random-number streams: one seed
// from the chain's generator, jumped i times for draw i (chunk_streams), so
// the result is the same whether they are made in parallel (parallel = true:
// RcppParallel's threads) or not.
class RatioEstimatingMCMC : public MetropolisMCMC
{

public:

  RatioEstimatingMCMC();

  RatioEstimatingMCMC(size_t number_of_iterations_in,
                      SymmetricProposalKernel* proposal_in,
                      RatioStep* step_in);

  RatioEstimatingMCMC(MCMCTermination* termination_in,
                      SymmetricProposalKernel* proposal_in,
                      RatioStep* step_in);

  virtual ~RatioEstimatingMCMC();

  RatioEstimatingMCMC(const RatioEstimatingMCMC &another);

  void operator=(const RatioEstimatingMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  MetropolisMCMC* metropolis_mcmc_duplicate() const;

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

  // points = N: the averaged acceptance ratio (above); parallel: the N draws
  // in RcppParallel's threads (the same results either way)
  EstimatedRatio& estimated_ratio() { return this->ratio; }
  void link_kernel(const std::vector<const MCMC*> &cycle) { if (this->ratio.get_step() != NULL) this->ratio.get_step()->link_kernel(cycle); }
  // (U5a) not as part of an annealed move's kernel (it moves the parameters
  // of a normalising constant, with an estimated ratio)
  bool reversible_parts(std::string &why) const
  {
    why = "an exchange move changes parameters, with an estimated acceptance ratio";
    return false;
  }

protected:

  StandardMCMCOutput* initialise_mcmc_output() const;

  // the decision for a proposed particle (after the uniform is drawn): the
  // move's factors first (a proposal they rule out draws no field), then the
  // fields' terms. Sets the acceptance probability, for tuning and adaptation.
  bool decide(RandomNumberGenerator &rng,
              double log_u,
              const Particle &particle,
              Particle &proposed_particle,
              double &acceptance_probability,
              Particle* refreshed = NULL,
              bool* refreshed_made = NULL) const;

  void accept(Particle &proposed_particle) const;

  void make_copy(const RatioEstimatingMCMC &another);

  // the estimated term (the step, the number of draws, ...)
  EstimatedRatio ratio;

};
}

#endif
