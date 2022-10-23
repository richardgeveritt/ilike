#ifndef PARAMETERESTIMATOR_H
#define PARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

class ParameterEstimator
{

public:

  ParameterEstimator();
  virtual ~ParameterEstimator();

  ParameterEstimator(const ParameterEstimator &another);

  void operator=(const ParameterEstimator &another);

protected:

  void make_copy(const ParameterEstimator &another);

};

#endif
