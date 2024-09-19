#ifndef ITERATIONSMCMCTERMINATION_H
#define ITERATIONSMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

// Checks to see

namespace ilike
{
class IterationsMCMCTermination : public MCMCTermination
{
  
public:
  
  IterationsMCMCTermination();
  
  IterationsMCMCTermination(size_t number_of_iterations_in);
  
  virtual ~IterationsMCMCTermination();
  
  IterationsMCMCTermination(const IterationsMCMCTermination &another);
  
  void operator=(const IterationsMCMCTermination &another);
  MCMCTermination* duplicate() const;
  
  bool terminate();
  
  void set_parameters(StandardMCMCOutput* mcmc_output);
  
protected:
  
  // not stored here
  size_t* counter;
  
  size_t number_of_iterations;
  
  void make_copy(const IterationsMCMCTermination &another);
  
};
}

#endif
