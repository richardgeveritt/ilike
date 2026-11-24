#ifndef ADAPTIVEITERATIONSMCMCTERMINATION_H
#define ADAPTIVEITERATIONSMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include "mcmc_termination.h"

namespace ilike
{

// The number of MCMC iterations of each SMC step chosen from the particles'
// acceptance rate (South, Pettitt & Drovandi 2019; used by Kerama, Thorne &
// Everitt, arXiv:2211.02172, §3.3.3; dev/plan.md 3.21 (c3)): the SMC first
// moves every particle once, estimates the acceptance rate p from that pass,
// and then makes R - 1 more iterations, R = ceil(log c / log(1 - p)), so that
// each particle moves at least once with probability about 1 - c (R at most
// max). The number for the current pass is shared by this termination and
// its copies in the chains (the SMC sets it before each pass).
class AdaptiveIterationsMCMCTermination : public MCMCTermination
{
public:

  AdaptiveIterationsMCMCTermination();
  AdaptiveIterationsMCMCTermination(double c_in, size_t maximum_in);
  virtual ~AdaptiveIterationsMCMCTermination();
  AdaptiveIterationsMCMCTermination(const AdaptiveIterationsMCMCTermination &another);
  void operator=(const AdaptiveIterationsMCMCTermination &another);
  MCMCTermination* duplicate() const;

  bool terminate();
  void set_parameters(StandardMCMCOutput* mcmc_output);
  // (not known in advance: no flat moves, which need it)
  size_t get_known_number_of_iterations() const { return 0; }

  // the iterations of the next pass (shared with the copies)
  void set_iterations(size_t iterations) { *this->iterations = iterations; }
  // R for an acceptance rate p of the pilot pass
  size_t iterations_for(double acceptance_rate) const;

protected:

  void make_copy(const AdaptiveIterationsMCMCTermination &another);

  size_t* counter = NULL;
  std::shared_ptr<size_t> iterations;
  double c = 0.2;
  size_t maximum = 100;
};

}

#endif
