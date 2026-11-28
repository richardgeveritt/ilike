#ifndef SMCPATH_H
#define SMCPATH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <vector>
#include "annealed_path.h"

namespace ilike
{

class SMCMCMCMove;
class Factors;
class Index;
class MCMC;
class IndependentProposalKernel;
struct ParticleTuning;

// The general engine for annealed paths (dev/plan.md 3.22 U4, "The engines"):
// the package's own SMC with MCMC moves, run on the inner model of an
// unnormalised distribution written as a block (its base measure and its
// terms). The normaliser's parameters follow the SMC's last schedule
// variable t from 0 to 1 (Sequencer::set_path), so the path is
//   vartheta(t) = start + t (end - start)
// at the schedule's points, with the SMC's moves, resampling, waste-free
// chains and so on. A run starts from the given states, which are draws
// from the distribution at start (set_initial_particles), and its evidence
// estimate is an estimate of log Z(end)/Z(start): unbiased for the ratio
// (in the sense of the exponential) with a fixed schedule and fixed kernels.
// single() uses an SMC with one particle (annealed importance sampling).
//
// It also draws from the distribution at given parameters by MCMC
// (simulate = sweeps(m)): m iterations of the inner model's kernel.
class SMCPath : public AnnealedPath
{
public:
  SMCPath();
  // smc_in: the path engine (owned; NULL if the path is not needed);
  // draw_factors_in, draw_index_in, draw_mcmc_in: the inner model's factors,
  // their index and its MCMC kernel, for draws (owned; NULL if not needed)
  SMCPath(SMCMCMCMove* smc_in,
          Factors* draw_factors_in,
          Index* draw_index_in,
          MCMC* draw_mcmc_in,
          const std::vector<std::string> &variables_in,
          const std::vector<std::string> &normaliser_in);
  ~SMCPath();
  SMCPath(const SMCPath &another);
  AnnealedPath* duplicate() const;

  double single(RandomNumberGenerator &rng, Parameters &state,
                const Parameters &start, const Parameters &end) const;
  // (threshold is the inner SMC's own resampling setting, not this argument)
  double population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                    const Parameters &start, const Parameters &end,
                    double threshold = 1.0) const;

  // (stage A) one particle started from the whole state (its variables other
  // than names, which the path moves from start to end) and returned at the
  // end of the path (names at end): AIS MCMC's path
  // tuning: the caller's particle's (esjd in SMC), carried by the path's
  // particle during its moves (A2)
  double single_state(RandomNumberGenerator &rng, Parameters &state,
                      const std::vector<std::string> &names,
                      const Parameters &start, const Parameters &end,
                      const ParticleTuning* tuning = NULL) const;

  // (A2) the path's moves are these moves of the caller's cycle (an annealed
  // move's kernel), not its own
  void link_kernel(const std::vector<const MCMC*> &moves);

  // the field variables of state after `sweeps` iterations of the kernel at
  // the parameters `at`, from their values in state
  void draw(RandomNumberGenerator &rng, Parameters &state, const Parameters &at, size_t sweeps) const;

  // a draw of the field variables from the base measure, given the
  // parameters it may read
  Parameters draw_base(RandomNumberGenerator &rng, const Parameters &given) const;

  // (normalise) built with an annealing sequence instead of a path: an SMC
  // from the base measure to the distribution at `at`, the potential raised
  // to a power from 0 to 1; its evidence estimate is log Z(at) relative to
  // the base measure
  double log_evidence(RandomNumberGenerator &rng, const Parameters &at) const;

  void set_base(IndependentProposalKernel* base_in);
  // the data the inner estimators were made with (kept alive here)
  void keep_data(std::shared_ptr<void> kept_in) { this->kept.push_back(kept_in); }

  bool has_path() const { return this->smc != NULL; }
  bool can_draw() const { return this->draw_mcmc != NULL; }

protected:
  double run_from(RandomNumberGenerator &rng, std::vector<Parameters> &points,
                  const std::vector<std::string> &names,
                  const Parameters &start, const Parameters &end,
                  const ParticleTuning* tuning = NULL) const;

  // (A2) A run writes to its SMC (its path, starting particles, generator
  // and sequencer), so runs on several threads at once each take their own
  // copy: a pool, the original first, with a copy made only when every
  // copy is in use (one per thread at most). Each run sets everything a
  // run reads, so which copy it gets does not change its result.
  SMCMCMCMove* acquire() const;
  void release(SMCMCMCMove* smc_in) const;
  struct Lease
  {
    const SMCPath* path;
    SMCMCMCMove* smc;
    Lease(const SMCPath* path_in) : path(path_in), smc(path_in->acquire()) {}
    ~Lease() { this->path->release(this->smc); }
  };
  mutable std::mutex pool_mutex;
  mutable std::condition_variable pool_freed;
  mutable std::vector<std::unique_ptr<SMCMCMCMove>> pool_copies; // (besides the original)
  mutable std::vector<SMCMCMCMove*> pool_idle;

  mutable std::unique_ptr<SMCMCMCMove> smc;
  mutable std::unique_ptr<Factors> draw_factors;
  std::unique_ptr<Index> draw_index;
  std::unique_ptr<MCMC> draw_mcmc;
  std::unique_ptr<IndependentProposalKernel> base;
  std::vector<std::string> variables;
  std::vector<std::string> normaliser;
  std::vector<std::shared_ptr<void>> kept;
};

// A pointer that copies its object (T::duplicate) when it is copied, so that
// objects copied for threads each have their own engine.
template<class T>
class ClonePtr
{
public:
  ClonePtr() {}
  explicit ClonePtr(T* p) : pointer(p) {}
  ClonePtr(const ClonePtr &another) : pointer(another.pointer ? static_cast<T*>(another.pointer->duplicate()) : NULL) {}
  ClonePtr& operator=(const ClonePtr &another)
  {
    if (this != &another)
      this->pointer.reset(another.pointer ? static_cast<T*>(another.pointer->duplicate()) : NULL);
    return *this;
  }
  T* operator->() const { return this->pointer.get(); }
  T* get() const { return this->pointer.get(); }
  explicit operator bool() const { return this->pointer != NULL; }
  void reset(T* p) { this->pointer.reset(p); }
private:
  std::unique_ptr<T> pointer;
};

}

#endif
