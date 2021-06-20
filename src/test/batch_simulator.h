#include <RcppArmadillo.h>
using namespace Rcpp;

#include "particles.h"

#ifndef BATCHSIMULATOR_H
#define BATCHSIMULATOR_H

class BatchSimulator
{
public:

  BatchSimulator(void);
  BatchSimulator(const BatchSimulator &another);
  virtual ~BatchSimulator(void);

  void operator=(const BatchSimulator &another);
  virtual BatchSimulator* duplicate() const=0;

  virtual Particles simulate(void) const=0;

protected:

  void make_copy(const BatchSimulator &another);

  // Add in function to do simuulation of parts that cannot be done in parallel.

};

#endif
