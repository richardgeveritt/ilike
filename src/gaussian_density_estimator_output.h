#ifndef GAUSSIANDENSITYESTIMATOROUTPUT_H
#define GAUSSIANDENSITYESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "density_estimator_output.h"

class VectorParameterEstimator;
class MatrixParameterEstimator;
class GaussianDensityEstimator;

class GaussianDensityEstimatorOutput : public DensityEstimatorOutput
{

public:

  GaussianDensityEstimatorOutput();
  
  GaussianDensityEstimatorOutput(GaussianDensityEstimator* estimator_in);
  
  virtual ~GaussianDensityEstimatorOutput();

  GaussianDensityEstimatorOutput(const GaussianDensityEstimatorOutput &another);

  void operator=(const GaussianDensityEstimatorOutput &another);
  DensityEstimatorOutput* duplicate() const;

  void fit(const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  double evaluate(const Data &point) const;
  
protected:
  
  // not stored heres
  GaussianDensityEstimator* estimator;
  
  // Stored here.
  VectorParameterEstimator* mean_estimator;
  MatrixParameterEstimator* covariance_estimator;

  void make_copy(const GaussianDensityEstimatorOutput &another);

};

#endif
