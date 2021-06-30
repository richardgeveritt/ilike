#ifndef SMCWORKER_H
#define SMCWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "particles.h"
#include "simulation.h"

class SMCWorker
{
public:

  SMCWorker(void);
  SMCWorker(const SMCWorker &another);
  virtual ~SMCWorker(void);

  void operator=(const SMCWorker &another);
  virtual SMCWorker* duplicate() const=0;

  virtual Particles simulate(void) const=0;

protected:

  void make_copy(const SMCWorker &another);

  // Add in function to do simuulation of parts that cannot be done in parallel.

};

#endif
