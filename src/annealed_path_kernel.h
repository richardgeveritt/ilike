#ifndef ANNEALEDPATHKERNEL_H
#define ANNEALEDPATHKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <string>
#include <vector>
#include "smc_path.h"
#include "index.h"

namespace ilike
{

class MCMC;
class Particle;

// What a state step with an annealed path keeps (dev/plan.md 3.22, "Every
// annealed path on the general engine"; AIS MCMC, C317-C318, and annealed
// reversible jump, stage B): the number of bridges, the path (the general
// engine, NULL with no bridges), the index of every factor (after a path,
// which changes the latent or jump variables, the new state is evaluated on
// all of them) and the cycle positions of the kernel's moves (the path runs
// the cycle's own moves, linked when the cycle is built or copied).
struct AnnealedPathKernel
{
  size_t bridges = 0;
  ClonePtr<SMCPath> path;
  // (stage B) annealed reversible jump: the path for a jump to the smaller
  // model, run from the larger model's side (the current state) with the
  // powers swapped; path is then the one for a jump to the larger model
  ClonePtr<SMCPath> reverse_path;
  std::unique_ptr<Index> state_index;
  std::vector<size_t> kernel;

  AnnealedPathKernel() {}
  AnnealedPathKernel(size_t bridges_in, SMCPath* path_in, Index* state_index_in, const std::vector<size_t> &kernel_in);
  AnnealedPathKernel(const AnnealedPathKernel &another);
  AnnealedPathKernel& operator=(const AnnealedPathKernel &another);

  bool runs() const { return (this->bridges > 0) && this->path; }
  // the kernel's moves from the cycle's (MCMC::link_kernel)
  void link_kernel(const std::vector<const MCMC*> &cycle, const std::string &function) const;
  // (A2) refused: an estimated factor in index (each bridge would draw a new
  // estimate); for state-space models this will be AIS with conditional SMC (U8)
  void refuse_estimates(const Particle &current, const Index* index, const std::string &function) const;
};

}

#endif
