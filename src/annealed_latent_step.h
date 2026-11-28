#ifndef ANNEALEDLATENTSTEP_H
#define ANNEALEDLATENTSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <vector>
#include "ratio_step.h"
#include "annealed_path_kernel.h"

namespace ilike
{

// AIS MCMC as a state step of the ratio-estimating moves (dev/plan.md 3.22
// U5a, C310; the common engine, C311; Karagiannis & Andrieu 2013; Yildirim,
// Andrieu & Doucet, arXiv 1809.02527, "MCMC AIS"). Given the current (theta,
// x) and the proposed theta': on the path vartheta_k = theta + (k/(K+1))(theta' - theta),
// k = 0, ..., K+1, x_0 = x and, for k = 1, ..., K, x_k from the latent kernel
// at vartheta_k (sweeps times a symmetric scan of the kernel's moves, a
// reversible kernel leaving pi(x | vartheta_k) invariant); the new state is
// (theta', x_K), and the target part of the log acceptance ratio is
//   sum_{k=0}^K log pi(vartheta_{k+1}, x_k) - log pi(vartheta_k, x_k).
// Valid for any K >= 0 (K = 0: theta moves with x fixed).
//
// The path is a run of the general engine ("Every annealed path on the
// general engine", stage A): the package's SMC with one particle (annealed
// importance sampling) on the model, along a parameter path in theta
// (Sequencer::set_path), its moves the kernel's (made from the model's move
// statements, as a symmetric scan), its evidence estimate the sum above. K =
// 0 needs no run: the factors are evaluated at theta'.
class AnnealedLatentStep : public RatioStep
{
public:

  // state_index_in: every factor (owned; the path changes the latent
  // variables, so the new state is evaluated on all of them)
  // kernel_in: the kernel's moves' positions in the cycle (A2: the path runs
  // the cycle's own moves, linked by link_kernel)
  AnnealedLatentStep(size_t bridges_in, SMCPath* path_in, Index* state_index_in,
                     const std::vector<size_t> &kernel_in = std::vector<size_t>());
  AnnealedLatentStep(const AnnealedLatentStep &another);
  virtual ~AnnealedLatentStep();
  RatioStep* duplicate() const;

  // (never called: a state step's draw is draw_state)
  double draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const;
  bool changes_state() const { return true; }
  double draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const;
  void link_kernel(const std::vector<const MCMC*> &cycle);

protected:

  AnnealedPathKernel annealing;
};

}

#endif
