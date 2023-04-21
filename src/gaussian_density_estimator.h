#ifndef GAUSSIANDENSITYESTIMATOR_H
#define GAUSSIANDENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "density_estimator.h"

class VectorParameterEstimator;
class MatrixParameterEstimator;
class DensityEstimatorOutput;

class GaussianDensityEstimator : public DensityEstimator
{

public:

  GaussianDensityEstimator();
  
  GaussianDensityEstimator(const std::vector<std::string> &variables_in);
  
  GaussianDensityEstimator(const std::vector<std::string> &variables_in,
                           bool unbiased_in);

  virtual ~GaussianDensityEstimator();

  GaussianDensityEstimator(const GaussianDensityEstimator &another);

  void operator=(const GaussianDensityEstimator &another);
  DensityEstimator* duplicate() const;
  
  DensityEstimatorOutput* initialise();
  
  bool get_unbiased() const;

  /*
  void fit(const std::vector<Parameters> &points,
           arma::colvec normalised_log_weights);
  
  double evaluate(const Data &point) const;
  */
  
protected:
  
  // Stored here.
  //VectorParameterEstimator* mean_estimator;
  //MatrixParameterEstimator* covariance_estimator;
  
  bool unbiased;

  void make_copy(const GaussianDensityEstimator &another);

};

#endif
