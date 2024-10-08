#ifndef SEMCMCTERMINATION_H
#define SEMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

namespace ilike
{
class SEMCMCTermination : public MCMCTermination
{
  
public:
  
  SEMCMCTermination();
  
  SEMCMCTermination(double threshold_in,
                    size_t max_number_of_iterations_in);
  
  virtual ~SEMCMCTermination();
  
  SEMCMCTermination(const SEMCMCTermination &another);
  
  void operator=(const SEMCMCTermination &another);
  MCMCTermination* duplicate() const;
  
  bool terminate();
  
  void set_parameters(StandardMCMCOutput* mcmc_output);
  
protected:
  
  // not stored here
  size_t* counter;
  
  double threshold;
  
  size_t max_number_of_iterations;
  
  void make_copy(const SEMCMCTermination &another);
  
};
}

#endif
