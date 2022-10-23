#ifndef SMCTERMINATION_H
#define SMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"

class SMCTermination
{

public:

  SMCTermination();
  virtual ~SMCTermination();

  SMCTermination(const SMCTermination &another);

  void operator=(const SMCTermination &another);
  virtual SMCTermination* duplicate() const=0;

  virtual bool terminate(double score)=0;

protected:

  void make_copy(const SMCTermination &another);

};

#endif
