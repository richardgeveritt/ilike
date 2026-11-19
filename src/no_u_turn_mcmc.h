#ifndef NOUTURNMCMC_H
#define NOUTURNMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <vector>

#include "metropolis_hastings_mcmc.h"

namespace ilike
{

class HMCProposalKernel;

// The No-U-Turn sampler (Hoffman & Gelman, arXiv:1111.4246v1, Algorithm 3,
// "efficient NUTS"; dev/plan.md 3.18 stage H3, C210). From the current
// point x0 with a fresh momentum r0 ~ N(0, C^-1) and a slice variable
// u ~ U(0, exp(L(x0) - K(r0))), the trajectory is doubled, each time in a
// random direction, until a sub-trajectory makes a U-turn,
//   (x+ - x-) . C r- < 0  or  (x+ - x-) . C r+ < 0
// (the momentum entering through the velocity C r, for a mass matrix
// M = C^-1), or the error in the Hamiltonian exceeds Delta_max = 1000 (a
// divergence), or the tree reaches the maximum depth. The next point is drawn
// among the states inside the slice, by the progressive sampling between
// subtrees of the paper (a new subtree's candidate replaces the current one
// with probability min(1, n'/n) at the top level, n''/(n' + n'') inside). L is
// the move's log density (its factors in the move's index, plus log|dx/dz|
// on a transformed scale); the Hamiltonian system (momenta, leapfrog steps,
// gradients, transform) is the HMCProposalKernel's, whose step size eps
// (eps^2 its scale) and mass (C) the adaptors adapt, with the statistic of
// Algorithm 6: the average of min(1, exp(H0 - H)) over the final doubling's
// states.
class NoUTurnMCMC : public MetropolisHastingsMCMC
{
public:

  NoUTurnMCMC();
  NoUTurnMCMC(size_t number_of_iterations_in,
              HMCProposalKernel* proposal_in);
  NoUTurnMCMC(MCMCTermination* termination_in,
              HMCProposalKernel* proposal_in);
  virtual ~NoUTurnMCMC();

  NoUTurnMCMC(const NoUTurnMCMC &another);
  void operator=(const NoUTurnMCMC &another);
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

  // one NUTS transition, with the adaptation statistic, the depth reached and
  // whether a divergence stopped it
  Particle transition(RandomNumberGenerator &rng,
                      const Particle &particle,
                      double &statistic,
                      size_t &depth,
                      bool &divergent) const;

protected:

  struct State
  {
    arma::colvec x;
    arma::colvec r;
    arma::colvec gradient;
    double log_density;
    std::shared_ptr<Particle> particle; // (NULL: the current particle)
  };

  struct Tree
  {
    State minus;
    State plus;
    State candidate;
    double n;
    bool s;
    double alpha;
    double n_alpha;
    bool divergent;
  };

  Tree build_tree(RandomNumberGenerator &rng,
                  const Particle &particle,
                  const Parameters &base,
                  const State &from,
                  double log_u,
                  int direction,
                  size_t depth,
                  double eps,
                  double joint0) const;
  bool no_u_turn(const State &minus, const State &plus) const;
  const HMCProposalKernel* hamiltonian() const;
};

}

#endif
