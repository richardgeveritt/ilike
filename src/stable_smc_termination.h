#ifndef STABLESMCTERMINATION_H
#define STABLESMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_termination.h"

// Checks to see

namespace ilike
{
class StableSMCTermination : public SMCTermination
{
  
public:
  
  StableSMCTermination();
  
  StableSMCTermination(size_t number_in_a_row_in,
                       double threshold_in);
  
  virtual ~StableSMCTermination();
  
  StableSMCTermination(const StableSMCTermination &another);
  
  void operator=(const StableSMCTermination &another);
  SMCTermination* duplicate() const;
  
  bool terminate(double score);
  
protected:
  
  size_t counter;
  
  size_t number_in_a_row;
  double threshold;
  
  void make_copy(const StableSMCTermination &another);
  
};
}

#endif
