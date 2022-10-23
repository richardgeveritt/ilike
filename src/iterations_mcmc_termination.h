#ifndef ITERATIONSMCMCTERMINATION_H
#define ITERATIONSMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

// Checks to see

class IterationsMCMCTermination : public MCMCTermination
{

public:

  IterationsMCMCTermination();
  
  IterationsMCMCTermination(size_t number_of_iterations_in,
                            size_t* counter_pointer);

  virtual ~IterationsMCMCTermination();

  IterationsMCMCTermination(const IterationsMCMCTermination &another);

  void operator=(const IterationsMCMCTermination &another);
  MCMCTermination* duplicate() const;

  bool terminate();
  
protected:
  
  // not stored here
  size_t* counter;
  
  size_t number_of_iterations;
  
  void make_copy(const IterationsMCMCTermination &another);

};

#endif
