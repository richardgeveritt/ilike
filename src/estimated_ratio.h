#ifndef ESTIMATEDRATIO_H
#define ESTIMATEDRATIO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include "ratio_step.h"

namespace ilike
{

// The estimated term of an acceptance ratio (dev/plan.md 3.22 U2, U1b),
// shared by the moves that estimate their ratios (RatioEstimatingMCMC with a
// symmetric proposal, RatioEstimatingMHMCMC with any proposal). One draw of
// the step, or N:
//  averaged (MHAAR, Andrieu, Yildirim, Doucet & Chopin, arXiv 2101.01253
//  Algorithm 1; exact): a direction c, equally likely 1 or 2;
//    c = 1: t_i = step.draw(theta, vartheta), i = 1..N; log term = log (1/N) sum_i exp(t_i);
//    c = 2: t_1 = -step.draw(theta, vartheta) (the reverse move's term for the
//           same u, since r o phi = 1/r), t_i = step.draw(vartheta, theta), i = 2..N;
//           log term = -log (1/N) sum_i exp(t_i);
//  noisy (Alquier, Friel, Everitt & Boland 2016, the noisy exchange; not
//  exact, its error falls as N grows): always the c = 1 average.
// For N = 1 both directions give the same term, so none is drawn. The N draws
// use their own streams (one seed from the chain's generator, jumped i times
// for draw i), so the result does not depend on whether they are made in
// parallel.
class EstimatedRatio
{
public:

  EstimatedRatio();
  EstimatedRatio(RatioStep* step_in);
  EstimatedRatio(const EstimatedRatio &another);
  EstimatedRatio& operator=(const EstimatedRatio &another);

  void set_points(size_t points_in) { this->points = std::max(points_in, size_t(1)); }
  void set_noisy(bool noisy_in) { this->noisy = noisy_in; }
  void set_parallel(bool parallel_in, size_t number_of_threads_in);

  double log_term(RandomNumberGenerator &rng,
                  const Parameters &current,
                  const Parameters &proposed) const;

protected:

  // owned
  std::unique_ptr<RatioStep> step;
  size_t points = 1;
  bool noisy = false;
  bool parallel = false;
  size_t number_of_threads = 1;

};

}

#endif
