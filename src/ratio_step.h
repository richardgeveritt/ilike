#ifndef RATIOSTEP_H
#define RATIOSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{

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

  // Whether a draw also gives new values of variables of the state (phi's
  // first part, a latent z' as in MHAAR-RB or reversible jump). Not written
  // yet: no step does so, and RatioEstimatingMCMC refuses one that does.
  virtual bool changes_state() const { return false; }

};

}

#endif
