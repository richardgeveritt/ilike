#ifndef DENSITYESTIMATOROUTPUT_H
#define DENSITYESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

//class SMCOutput;

class DensityEstimatorOutput
{

public:

  DensityEstimatorOutput();
  virtual ~DensityEstimatorOutput();

  DensityEstimatorOutput(const DensityEstimatorOutput &another);

  void operator=(const DensityEstimatorOutput &another);
  virtual DensityEstimatorOutput* duplicate() const=0;

  virtual void fit(const std::vector<Parameters> &points,
                   const arma::colvec &normalised_log_weights)=0;
  void fit(const std::vector<Parameters> &points);
  
  virtual double evaluate(const Parameters &point) const=0;

protected:
  
  size_t n;

  void make_copy(const DensityEstimatorOutput &another);

};

#endif
