#ifndef EARLYREJECTION_H
#define EARLYREJECTION_H

#include "particle.h"
#include "index.h"

namespace ilike
{

// The Metropolis-Hastings decision with early rejection (dev/plan.md 3.21),
// for a proposed particle with factors whose draws are deferred: true if the
// proposal is accepted. log_u is the uniform's log, drawn before the deferred
// draws; extra the terms of the log ratio other than numerator - denominator
// (the proposal kernels, a transform's Jacobian). The proposal is accepted iff
// log_u < numerator - denominator + extra, i.e. iff the numerator exceeds
// threshold = log_u + denominator - extra; the numerator is evaluated only as
// far as needed to tell (so the decision is the one a full evaluation makes,
// and a proposal rejected early never makes its deferred draws).
inline bool early_rejection_accepts(double log_u,
                                    double denominator,
                                    double extra,
                                    Particle &proposed,
                                    const Index* index)
{
  if (denominator == -arma::datum::inf)
    return false;
  bool below = false;
  double numerator = proposed.evaluate_likelihoods_with_threshold(index, log_u + denominator - extra, below);
  if (below || (numerator == -arma::datum::inf))
    return false;
  return log_u < numerator - denominator + extra;
}

}

#endif
