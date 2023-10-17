#ifndef ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"

class EnKGaussianMeasurementCovarianceEstimator : public GaussianMeasurementCovarianceEstimator
{

public:

  EnKGaussianMeasurementCovarianceEstimator();

  virtual ~EnKGaussianMeasurementCovarianceEstimator();

  EnKGaussianMeasurementCovarianceEstimator(const EnKGaussianMeasurementCovarianceEstimator &another);

  void operator=(const EnKGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  arma::mat get_Cygivenx() const;
  
  arma::mat get_measurement_covariance() const;
  
  void setup();
  void setup(const Parameters &parameters);

protected:
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);

  void make_copy(const EnKGaussianMeasurementCovarianceEstimator &another);

};

#endif
