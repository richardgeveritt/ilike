#ifndef RATIOSTEP_H
#define RATIOSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"
#include <vector>

namespace ilike
{

class Particle;
class Index;
class MCMC;

// One auxiliary draw of a move that estimates its acceptance ratio
// (dev/plan.md 3.22 A, U2). In the involution view (2101.01253 §2.1; 2012.14881)
// the move draws variables u ~ Q(from -> to) beside the proposed parameters,
// and its acceptance ratio is r_u(from, to) = (the ratio of the move's own
// factors) x (a term depending on u). draw makes one such u and returns the
// log of that term. Since the involution gives r o phi = 1/r, the term of the
// reverse move with the same u is minus the returned value, which the
// averaged acceptance ratio uses (RatioEstimatingMCMC). A step must make all
// its random numbers from the generator it is given, and must be safe to call
// from several threads at once (it may read, not write, shared state).
class RatioStep
{
public:

  virtual ~RatioStep() {}

  virtual RatioStep* duplicate() const=0;

  // u ~ Q(from -> to); the log of r_u(from, to)'s u-dependent term
  // (-infinity when u cannot be drawn there, e.g. outside the support)
  virtual double draw(RandomNumberGenerator &rng,
                      const Parameters &from,
                      const Parameters &to) const=0;

  // (dev/plan.md 3.22, the common engine, C311) A state step: one whose draw
  // also gives new values of variables of the state (phi's first part: the
  // latent x_K of AIS MCMC, a latent z' of MHAAR-RB, reversible jump).
  // draw_state is given the current particle and the proposal's particle
  // (the proposed parameters, the other variables as they were); it draws the
  // step's auxiliary variables, makes `proposed` the new state, with the
  // factors in `index` (the move's) evaluated, and returns the log of the
  // target part of the acceptance ratio (the move adds the proposal's own
  // terms). -infinity when the new state is outside the support.
  virtual bool changes_state() const { return false; }
  virtual double draw_state(RandomNumberGenerator &rng,
                            const Index* index,
                            const Particle &current,
                            Particle &proposed) const
  {
    Rcpp::stop("RatioStep::draw_state: this step does not change the state.");
  }
  // (C312) A state step that uses the direction c of the averaged
  // construction (MHAAR-RB: c = 1 or 2, drawn by the move with probability
  // 1/2 each), and one that, besides the proposal, makes a refreshed current
  // state for a rejection (MHAAR-RB's optional stage; *refreshed_made true
  // when it did). By default: the plain draw_state.
  virtual bool uses_direction() const { return false; }
  virtual bool refreshes() const { return false; }
  virtual double draw_state(RandomNumberGenerator &rng,
                            const Index* index,
                            const Particle &current,
                            Particle &proposed,
                            int direction,
                            Particle* refreshed,
                            bool* refreshed_made) const
  {
    return this->draw_state(rng, index, current, proposed);
  }
  // (C313) A state step that can be averaged over points = N draws with a
  // selection among the N new states (RmJ-MCMC, 2101.01253 Algorithm 3), and
  // the probability omega(from -> to, 1) of the direction c = 1 for a move
  // from `from` to `to` (Remark "choosing Q1 and Q2 with different
  // probabilities"; 1/2 by default)
  virtual bool averages() const { return false; }
  virtual double first_direction_probability(const Particle &from, const Particle &to) const { return 0.5; }
  // (A2) a step that runs moves of the cycle (AIS MCMC's kernel): see MCMC::link_kernel
  virtual void link_kernel(const std::vector<const MCMC*> &cycle) {}

};

}

#endif
