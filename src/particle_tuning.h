#ifndef PARTICLETUNING_H
#define PARTICLETUNING_H

#include <cstddef>
#include <utility>
#include <vector>

namespace ilike
{

class ESJDAdaptor;
class Particle;
class ProposalKernel;

// One particle's own tuning values in SMC (esjd, dev/plan.md 3.19, C211):
// for each tuned kernel (by its instance index), a pointer to the particle's
// column of values in the kernel's ESJDAdaptor, and the particle's score so
// far in this move step. Owned by the SMC move, one per particle; a particle
// points at it (Particle::tuning) only while it is moved, and only the thread
// moving that particle touches it.
struct ParticleTuning
{
  struct Row
  {
    int instance;
    const double* values;
    const ESJDAdaptor* adaptor;
    // the Rao-Blackwellised expected squared jumps alpha |x~ - x|^2_M summed
    // over the chain's Metropolis(-Hastings) iterations, and their number (C213)
    mutable double jumps = 0.0;
    mutable size_t iterations = 0;
    // choose's row (C224): every proposal of the move counts towards it
    bool choice = false;
  };
  std::vector<Row> rows;

  bool has_choice() const
  {
    for (const auto &row : this->rows)
      if (row.choice)
        return true;
    return false;
  }

  const double* for_kernel(int instance) const
  {
    for (const auto &row : this->rows)
      if (row.instance == instance)
        return row.values;
    return 0;
  }
};

// After a Metropolis(-Hastings) proposal from current to proposed by
// `proposal` (or a kernel inside it) with acceptance probability alpha: add
// alpha |x~ - x|^2_M to the rows of the kernels tuned within `proposal`
// (Fearnhead & Taylor's score; esjd_adaptor.cpp). Call only when
// current.tuning is set.
void observe_tuned_proposal(const Particle &current,
                            const Particle &proposed,
                            const ProposalKernel* proposal,
                            double alpha);

// alpha = min(1, exp(log_ratio)), 0 for NaN
inline double acceptance_from_log_ratio(double log_ratio)
{
  if (log_ratio >= 0.0)
    return 1.0;
  if (log_ratio != log_ratio)
    return 0.0;
  return __builtin_exp(log_ratio);
}

}

#endif
