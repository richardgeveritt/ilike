#include "enk_gaussian_measurement_covariance_estimator_output.h"

EnKGaussianMeasurementCovarianceEstimatorOutput::EnKGaussianMeasurementCovarianceEstimatorOutput()
  :GaussianMeasurementCovarianceEstimatorOutput()
{
}

EnKGaussianMeasurementCovarianceEstimatorOutput::EnKGaussianMeasurementCovarianceEstimatorOutput(EnKGaussianMeasurementCovarianceEstimator* enk_estimator_in)
:GaussianMeasurementCovarianceEstimatorOutput()
{
  this->enk_estimator = enk_estimator_in;
}

EnKGaussianMeasurementCovarianceEstimatorOutput::~EnKGaussianMeasurementCovarianceEstimatorOutput()
{
}

EnKGaussianMeasurementCovarianceEstimatorOutput::EnKGaussianMeasurementCovarianceEstimatorOutput(const EnKGaussianMeasurementCovarianceEstimatorOutput &another)
  :GaussianMeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void EnKGaussianMeasurementCovarianceEstimatorOutput::operator=(const EnKGaussianMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;

  GaussianMeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimatorOutput::duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimatorOutput(*this));
}

GaussianMeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimatorOutput::gaussian_duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimatorOutput(*this));
}

void EnKGaussianMeasurementCovarianceEstimatorOutput::make_copy(const EnKGaussianMeasurementCovarianceEstimatorOutput &another)
{
}

/*
arma::rowvec EnKGaussianMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance()
{
  return arma::rowvec();
}

arma::rowvec EnKGaussianMeasurementCovarianceEstimatorOutput::get_measurement_random_shift()
{
  return arma::rowvec();
}
*/

void EnKGaussianMeasurementCovarianceEstimatorOutput::simulate(const Parameters &parameters)
{
  
}

arma::mat EnKGaussianMeasurementCovarianceEstimatorOutput::get_measurement_covariance()
{
  return arma::mat();
}
