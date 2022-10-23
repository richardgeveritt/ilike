#ifndef GAUSSIANDENSITYESTIMATOR_H
#define GAUSSIANDENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "density_estimator.h"

class VectorParameterEstimator;
class MatrixParameterEstimator;

class GaussianDensityEstimator : public DensityEstimator
{

public:

  GaussianDensityEstimator();
  
  GaussianDensityEstimator(bool unbiased_in);

  virtual ~GaussianDensityEstimator();

  GaussianDensityEstimator(const GaussianDensityEstimator &another);

  void operator=(const GaussianDensityEstimator &another);
  DensityEstimator* duplicate() const;

  void fit(const std::vector<Parameters> &points,
           arma::colvec normalised_log_weights);
  
  double evaluate(const Parameters &point) const;
  
protected:
  
  // Stored here.
  VectorParameterEstimator* mean_estimator;
  MatrixParameterEstimator* covariance_estimator;
  
  bool unbiased;

  void make_copy(const GaussianDensityEstimator &another);

};

#endif
