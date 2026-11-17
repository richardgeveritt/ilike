#ifndef MCMCTERMINATION_H
#define MCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"
#include "standard_mcmc_output.h"

namespace ilike
{

class MCMC;

class MCMCTermination
{
  
public:
  
  MCMCTermination();
  virtual ~MCMCTermination();
  
  MCMCTermination(const MCMCTermination &another);
  
  void operator=(const MCMCTermination &another);
  virtual MCMCTermination* duplicate() const=0;
  
  virtual bool terminate()=0;
  
  virtual void set_parameters(StandardMCMCOutput* mcmc_output)=0;

  // The number of iterations a run will make, if known in advance (0 if
  // not): used to reserve the storage for the chain.
  virtual size_t get_known_number_of_iterations() const { return 0; }
  
protected:
  
  void make_copy(const MCMCTermination &another);
  
};
}

#endif
