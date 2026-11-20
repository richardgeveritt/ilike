#ifndef DUALAVERAGINGADAPTOR_H
#define DUALAVERAGINGADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>

#include "mcmc_adaptor.h"

namespace ilike
{

class HMCProposalKernel;

// The step size of hmc or nuts by dual averaging (Hoffman & Gelman,
// arXiv:1111.4246v1, Algorithms 5 and 6; Nesterov 2009; dev/plan.md 3.18
// stage H4, C219). With the acceptance statistic alpha_m of iteration m
// (hmc's acceptance probability; nuts's average over the last doubling):
//   H_m = (1 - 1/(m + t0)) H_{m-1} + (delta - alpha_m)/(m + t0)
//   log eps_m = mu - sqrt(m)/gamma H_m
//   log epsbar_m = m^-kappa log eps_m + (1 - m^-kappa) log epsbar_{m-1}
// with mu = log(10 eps_0), eps_0 from FindReasonableEpsilon (Algorithm 4:
// from the chain's starting point and one momentum, the step doubled or
// halved until the acceptance probability of one leapfrog step crosses 1/2)
// at the first adaptation step, or the step given (find_step = false). The
// iterates eps_m are not meant to converge; the average epsbar_m is. So it
// is used only with until = criterion: when the criterion holds, the step is
// fixed at epsbar (on_stop), and the chain from then on is a homogeneous
// Markov chain. The kernel's step is eps, its scale eps^2.
// With mass = diagonal or dense (stage H5, C220), the inverse mass matrix C
// is also estimated, in Stan's windows over a warm-up of W iterations (W from
// until = iterations(W)): a first fast window (75 iterations for W = 1000)
// adapting the step only; slow windows of 25, 50, 100, ... iterations (the
// last stretched to the end of the slow part), each estimating the
// covariance of the chain on the move's scale (regularised towards a small
// multiple of the identity, as Stan: n/(n + 5) S + 1e-3 5/(n + 5) I) and
// setting C to it, after which the step is found again
// (FindReasonableEpsilon) and dual averaging restarted; a last fast window
// (50) adapting the step only.
class DualAveragingAdaptor : public MCMCAdaptor
{

public:

  DualAveragingAdaptor();
  // mass: "none" (the mass given is kept), "diagonal" or "dense" (C220)
  DualAveragingAdaptor(double target_in, double gamma_in, double t0_in, double kappa_in, bool find_step_in,
                       const std::string &mass_in = "none");
  virtual ~DualAveragingAdaptor();

  DualAveragingAdaptor(const DualAveragingAdaptor &another);
  void operator=(const DualAveragingAdaptor &another);
  MCMCAdaptor* duplicate() const;

  void adapt(const AdaptationStep &step);
  void on_stop();

  void append_state(std::vector<std::string> *names,
                    std::vector<double> *values,
                    bool final,
                    const std::string &prefix) const;

  bool watched(const std::string &name, std::vector<double> &values) const;
  std::vector<std::string> watchable() const;
  bool acceptance_target(double &target) const;
  std::string state_label() const;

  // Algorithm 4 from the particle's point, starting from the step eps
  static double find_reasonable_epsilon(RandomNumberGenerator &rng,
                                        const HMCProposalKernel &kernel,
                                        const Particle &particle,
                                        double eps);

protected:

  void make_copy(const DualAveragingAdaptor &another);
  HMCProposalKernel* kernel() const;
  void set_step(double log_step);

  double target;
  double gamma;
  double t0;
  double kappa;
  bool find_step;
  std::string mass;

  // the windows (mass): the iterations ending each slow window
  std::vector<size_t> window_ends;
  size_t slow_start = 0;  // iterations before the first slow window
  size_t iteration = 0;   // adaptation steps so far
  size_t window_n = 0;    // points in the current window
  arma::colvec window_mean;
  arma::mat window_scatter;
  void plan_windows(size_t W);
  void restart(RandomNumberGenerator* rng, const Particle* particle);

  bool started = false;
  size_t m = 0;
  double mu = 0.0;
  double h_bar = 0.0;
  double log_eps = 0.0;
  double log_eps_bar = 0.0;
};

}

#endif
