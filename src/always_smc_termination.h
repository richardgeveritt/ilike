#ifndef ALWAYSSMCTERMINATION_H
#define ALWAYSSMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_termination.h"

namespace ilike
{
class AlwaysSMCTermination : public SMCTermination
{
  
public:
  
  AlwaysSMCTermination();
  
  virtual ~AlwaysSMCTermination();
  
  AlwaysSMCTermination(const AlwaysSMCTermination &another);
  
  void operator=(const AlwaysSMCTermination &another);
  SMCTermination* duplicate() const;
  
  bool terminate(double score);
  
protected:
  
  void make_copy(const AlwaysSMCTermination &another);
  
};
}

#endif
