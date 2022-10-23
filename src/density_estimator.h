#ifndef DENSITYESTIMATOR_H
#define DENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

//class SMCOutput;

class DensityEstimator
{

public:

  DensityEstimator();
  virtual ~DensityEstimator();

  DensityEstimator(const DensityEstimator &another);

  void operator=(const DensityEstimator &another);
  virtual DensityEstimator* duplicate() const=0;

  virtual void fit(const std::vector<Parameters> &points,
                   arma::colvec normalised_log_weights)=0;
  void fit(const std::vector<Parameters> &points);
  
  virtual double evaluate(const Parameters &point) const=0;

protected:
  
  size_t n;

  void make_copy(const DensityEstimator &another);

};

#endif
