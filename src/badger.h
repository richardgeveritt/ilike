#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

#ifndef BADGER_H
#define BADGER_H

class Badger
{
public:

  Badger(void);
  virtual ~Badger(void);

  Parameters thing;

protected:

};

#endif
