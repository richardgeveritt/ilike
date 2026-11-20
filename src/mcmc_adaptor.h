#ifndef MCMCADAPTOR_H
#define MCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"
#include "adaptation_criterion.h"

namespace ilike
{

class ProposalKernel;

// What one iteration of a Metropolis(-Hastings) move gives an adaptor
// (dev/plan.md 3.2, A1): the state before the move X_i, the proposed point
// Y_{i+1}, the state after it X_{i+1} (one of the two), and the acceptance
// probability alpha(X_i, Y_{i+1}) = min(1, ratio) (0 if either density is 0).
struct AdaptationStep
{
  const Particle *previous;
  const Particle *proposed;
  const Particle *current;
  double acceptance_probability;
  // the chain's generator (set by the hmc and nuts moves, for
  // dual_averaging's FindReasonableEpsilon, C219; NULL otherwise)
  RandomNumberGenerator* rng = NULL;
};

// Adapts the parameters of the proposal kernel that owns it, from the chain
// (Andrieu & Thoms 2008). Each chain adapts its own copy of the kernel (C98),
// so the adaptor is copied with the kernel and then pointed at the copy.
class MCMCAdaptor
{

public:

  MCMCAdaptor();
  virtual ~MCMCAdaptor();

  MCMCAdaptor(const MCMCAdaptor &another);
  void operator=(const MCMCAdaptor &another);
  virtual MCMCAdaptor* duplicate() const=0;

  // The kernel adapted: called when the adaptor is attached, and again on the
  // kernel's copy after the kernel is duplicated.
  virtual void set_proposal(ProposalKernel* proposal_in);

  // One step of adaptation after each iteration of the move.
  virtual void adapt(const AdaptationStep &step)=0;

  // The adaptation state written to the output (dev/plan.md 3.2): with
  // final false, the trace written at every stored iteration (scales and
  // means: O(d)); with final true, the state at the end of a chain (also
  // covariances and directions). Appends the column names (names != NULL)
  // or the values (values != NULL); prefix goes before each name.
  virtual void append_state(std::vector<std::string> *names,
                            std::vector<double> *values,
                            bool final,
                            const std::string &prefix) const=0;

  // Called before each move (some adaptors point the kernel at their own
  // parameters, which a copy of the kernel does not yet know about).
  virtual void prepare_for_move();

  // Moves that do not report an AdaptationStep cannot be adapted.
  void mcmc_adapt(const Particle &latest_particle,
                  size_t iteration_counter);

  // Stopping (dev/plan.md 3.18, item 10): adapt (one step) unless the
  // criterion given with until = has held; then check it.
  void adapt_unless_stopped(const AdaptationStep &step);
  void set_until(AdaptationCriterion* until_in); // (takes ownership)
  bool has_until() const { return this->until != NULL; }
  bool has_stopped() const { return this->stopped; }
  size_t stopped_at() const { return this->stopped_at_step; }
  // the number of adaptation steps before it stops, if until = iterations(n)
  // gives it (0 otherwise; C220)
  size_t warm_up_length() const { return (this->until != NULL) ? this->until->known_iterations() : 0; }

  // What a stopping criterion can read. watched(name, values): the current
  // values of a quantity the adaptor adapts ("scale", "covariance" (the
  // geometric mean of the eigenvalues, det^(1/d), and the mean variance,
  // trace/d), "covariance_elements", "mean", "component_scales"); FALSE if
  // it has no such quantity. default_watch(): the quantities it adapts.
  // acceptance_target(target): FALSE if it has no target rate.
  virtual bool watched(const std::string &name, std::vector<double> &values) const;
  virtual std::vector<std::string> watchable() const;
  virtual std::vector<std::string> default_watch() const;
  virtual bool acceptance_target(double &target) const;
  // how its output columns name what it adapts (for messages and the
  // column giving when it stopped)
  virtual std::string state_label() const;
  // called once, when the criterion given with until = holds (dual
  // averaging then fixes the step at its average, C219)
  virtual void on_stop();

protected:

  void make_copy(const MCMCAdaptor &another);

  // Not owned: the kernel that owns this adaptor.
  ProposalKernel* proposal;

  std::unique_ptr<AdaptationCriterion> until;
  bool stopped = false;
  size_t steps = 0;
  size_t stopped_at_step = 0;

};

}

#endif
