#ifndef ESJDADAPTOR_H
#define ESJDADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "smc_adaptor.h"
#include "particle_tuning.h"

namespace ilike
{

class Transform;

// Tuning of a kernel's parameters in SMC from how far the particles move
// (Fearnhead & Taylor 2013, arXiv:1005.1193v2; Buchholz, Chopin & Jacob 2021,
// arXiv:1808.07730v2, "FT"; dev/plan.md 3.19, C211). Each particle i moves
// with its own values h_i (a multiplier of a random walk's or Langevin's step;
// hmc's (step, steps); nuts's step), which the kernel reads from the particle
// (Particle::tuning). After the move, h_i is scored by the squared jump of
// its particle in the metric M^-1 (the inverse of the particles' variances,
// or covariance, on the kernel's scale), divided by the number of leapfrog
// steps for hmc (its cost): Lambda_i = |x_i' - x_i|^2_M / cost. Before the
// next move the h_i are resampled with weights a + Lambda_i (a >= 0; a = 0 by
// default, as Buchholz et al.), each perturbed (a continuous value multiplied
// by exp(perturb Z), Z ~ N(0, 1); a number of steps moved by -1, 0 or +1,
// each with probability 1/3, at least 1) and given to the particles in a
// random order. The first values are uniform on the ranges given. The jump is
// that of the particle's whole chain (with one MCMC iteration per SMC step,
// the jump of one step, as in the papers). The values in each column are one
// particle's (contiguous, for the particle's pointer).
class ESJDAdaptor : public SMCAdaptor
{
public:

  ESJDAdaptor();
  ESJDAdaptor(const std::string &variable_in,
              const std::vector<std::string> &names_in,
              const std::vector<double> &lower_in,
              const std::vector<double> &upper_in,
              const std::vector<bool> &integer_in,
              double constant_in,
              double perturb_in,
              bool full_metric_in,
              int cost_index_in);
  virtual ~ESJDAdaptor();

  ESJDAdaptor(const ESJDAdaptor &another);
  void operator=(const ESJDAdaptor &another);
  SMCAdaptor* duplicate() const;

  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);

  // before the SMC move's moves: the values for the N particles (new ones, or
  // resampled and perturbed from the last move's scores) and the metric
  virtual void update(RandomNumberGenerator &rng, SMCOutput* current_state, size_t number_of_particles);
  const double* column(size_t i) const { return this->values.colptr(i); }
  // |after - before|^2_M on the kernel's scale (read-only: called from the
  // threads moving the particles)
  double squared_jump(const Parameters &before, const Parameters &after) const;
  // after particle i's chain: the mean of alpha |x~ - x|^2_M over its
  // Metropolis(-Hastings) iterations, accumulated in its row of `tuning`
  // (C213); for a move with no such iterations (nuts), the realised jump
  // |after - before|^2_M of the chain; divided by the cost for hmc
  virtual void score(size_t i, const Parameters &before, const Parameters &after, const ParticleTuning &tuning);
  // the values' summary, for the output: per parameter, mean, sd, 10%, 50%, 90%
  virtual void summary(std::vector<std::string> &names, std::vector<double> &numbers) const;
  // pretune (C223): trial moves from each particle's starting state, after
  // update and before the move (start(i): particle i's state, a copy)
  virtual bool needs_trial_moves() const { return false; }
  // choose's kernel index (C224)
  virtual bool is_choice() const { return false; }
  virtual void trial_moves(RandomNumberGenerator &rng, size_t number_of_particles,
                           const std::function<Particle(size_t)> &start) {}
  // |jump|^2 in the metric (a jump on the kernel's scale)
  double mahalanobis(const arma::colvec &jump) const;
  // the metric's entry for element k of the block (a jump of one element, C226)
  double metric_at(size_t k) const
  {
    if ((this->metric.n_rows > 1) && (this->metric.n_cols > 1))
      return (k < this->metric.n_rows) ? this->metric(k, k) : 1.0;
    return (k < this->metric.n_elem) ? this->metric[k] : 1.0;
  }
  const std::string &get_variable() const { return this->variable; }
  size_t get_updates() const { return this->updates; }

protected:

  void make_copy(const ESJDAdaptor &another);
  void initialise(RandomNumberGenerator &rng, size_t number_of_particles);
  // value k of a resampled column, perturbed (C224: a hook for choose's index)
  virtual double perturbed(RandomNumberGenerator &rng, size_t k, double v) const;
  // the metric from the particles' variances (or covariance) on the kernel's scale
  void compute_metric(SMCOutput* current_state);

  std::string variable;
  std::vector<std::string> names;
  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<bool> integer;
  double constant;
  double perturb;
  bool full_metric;
  int cost_index;

  arma::mat values;       // (one column per particle)
  arma::colvec scores;    // (NaN: not scored)
  arma::mat metric;       // M^-1 (full) or its diagonal (a column)
  size_t updates = 0;     // (move steps so far, for the output)
};

// Several SMC adaptors on one proposal, applied in order
// (smc_adapt = [gaussian(2.38), esjd(…)], C211).
class CompositeSMCAdaptor : public SMCAdaptor
{
public:

  CompositeSMCAdaptor();
  CompositeSMCAdaptor(const std::vector<SMCAdaptor*> &parts_in); // (takes ownership)
  virtual ~CompositeSMCAdaptor();

  CompositeSMCAdaptor(const CompositeSMCAdaptor &another);
  void operator=(const CompositeSMCAdaptor &another);
  SMCAdaptor* duplicate() const;

  void set_proposal(ProposalKernel* proposal_in);
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  std::vector<SMCAdaptor*> get_parts() const;

protected:

  void make_copy(const CompositeSMCAdaptor &another);
  std::vector<std::unique_ptr<SMCAdaptor>> parts;
};

}

#endif
