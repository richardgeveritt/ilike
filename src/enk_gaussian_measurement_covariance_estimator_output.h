#ifndef ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator_output.h"

class EnKGaussianMeasurementCovarianceEstimator;

class EnKGaussianMeasurementCovarianceEstimatorOutput : public GaussianMeasurementCovarianceEstimatorOutput
{

public:

  EnKGaussianMeasurementCovarianceEstimatorOutput();
  
  EnKGaussianMeasurementCovarianceEstimatorOutput(EnKGaussianMeasurementCovarianceEstimator* enk_estimator_in);

  virtual ~EnKGaussianMeasurementCovarianceEstimatorOutput();

  EnKGaussianMeasurementCovarianceEstimatorOutput(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);

  void operator=(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const;
  
  //arma::rowvec get_measurement_state_for_covariance() const;
  //arma::rowvec get_measurement_random_shift();
  
  void specific_simulate(const Parameters &parameters);
  void subsample_specific_simulate(const Parameters &parameters);
  
  MeasurementCovarianceEstimator* get_estimator();
  GaussianMeasurementCovarianceEstimator* get_gaussian_estimator();

protected:
  
  // not stored here
  EnKGaussianMeasurementCovarianceEstimator* enk_estimator;

  void make_copy(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);

};

#endif
