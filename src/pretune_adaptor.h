#ifndef PRETUNEADAPTOR_H
#define PRETUNEADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "esjd_adaptor.h"

namespace ilike
{

// Pre-tuning of hmc in SMC (Buchholz, Chopin & Jacob 2021, arXiv:1808.07730v2,
// Algorithm "PR"; dev/plan.md 3.19 stage S3, C223). Before each move, every
// particle makes a trial HMC move from its starting state with its own
// eps_i ~ U(0, eps*] and L_i ~ U{1, ..., L_max}, which is then discarded,
// recording the energy error dE_i (the log acceptance ratio) and the score
//   Lambda_i = |x~_i - x_i|^2_M min(1, exp(dE_i)) / L_i
// (M as esjd's). Then:
//  - eps* from the median (least absolute deviations) regression
//    |dE| = a0 + a1 eps^2, solved for |log target| (target 0.9: the
//    acceptance rate about 90%); if the fit gives no solution, eps* is halved
//    (a0 already too large) or doubled (no growth with eps);
//  - each particle's (eps, L) for the real move drawn from the trial pairs
//    with weights Lambda (systematic resampling, random assignment);
//  - L_max raised by steps_change (5) when more than a share raise_share
//    (0.5) of the chosen L are above raise_above (0.8) L_max, lowered by
//    steps_change (to at least steps_change) when more than lower_share (0.9) are below
//    lower_below (0.5) L_max (the paper: "a large proportion close to / far
//    away from L_max"; these numbers are the defaults, C225).
// The trial moves target the current (reweighted) distribution, which the
// real move targets (the paper's trial flow targets the previous one).
class PretuneAdaptor : public ESJDAdaptor
{
public:

  PretuneAdaptor();
  PretuneAdaptor(const std::string &variable_in, double first_max_step_in, size_t first_max_steps_in,
                 double target_in, bool full_metric_in,
                 double raise_share_in = 0.5, double raise_above_in = 0.8,
                 double lower_share_in = 0.9, double lower_below_in = 0.5, size_t steps_change_in = 5);
  virtual ~PretuneAdaptor();

  PretuneAdaptor(const PretuneAdaptor &another);
  void operator=(const PretuneAdaptor &another);
  SMCAdaptor* duplicate() const;

  void update(RandomNumberGenerator &rng, SMCOutput* current_state, size_t number_of_particles);
  bool needs_trial_moves() const { return true; }
  void trial_moves(RandomNumberGenerator &rng, size_t number_of_particles,
                   const std::function<Particle(size_t)> &start);
  void score(size_t i, const Parameters &before, const Parameters &after, const ParticleTuning &tuning) {}
  void summary(std::vector<std::string> &names, std::vector<double> &numbers) const;

  // the median regression |dE| = a0 + a1 x (x = eps^2), by iteratively
  // reweighted least squares
  static void median_regression(const arma::colvec &x, const arma::colvec &y, double &a0, double &a1);

protected:

  void make_copy(const PretuneAdaptor &another);

  double max_step;   // eps*
  size_t max_steps;  // L_max
  double target;
  double raise_share;
  double raise_above;
  double lower_share;
  double lower_below;
  size_t steps_change;
};

}

#endif
