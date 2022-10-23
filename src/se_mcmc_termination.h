#ifndef SEMCMCTERMINATION_H
#define SEMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

class SEMCMCTermination : public MCMCTermination
{

public:

  SEMCMCTermination();
  
  SEMCMCTermination(double threshold_in);

  virtual ~SEMCMCTermination();

  SEMCMCTermination(const SEMCMCTermination &another);

  void operator=(const SEMCMCTermination &another);
  MCMCTermination* duplicate() const;

  bool terminate();
  
protected:
  
  double threshold;
  
  void make_copy(const SEMCMCTermination &another);

};

#endif
