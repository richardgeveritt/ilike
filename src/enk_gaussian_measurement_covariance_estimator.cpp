#include "enk_gaussian_measurement_covariance_estimator.h"
#include "enk_gaussian_measurement_covariance_estimator_output.h"

EnKGaussianMeasurementCovarianceEstimator::EnKGaussianMeasurementCovarianceEstimator()
  :GaussianMeasurementCovarianceEstimator()
{
}
EnKGaussianMeasurementCovarianceEstimator::~EnKGaussianMeasurementCovarianceEstimator()
{
}

EnKGaussianMeasurementCovarianceEstimator::EnKGaussianMeasurementCovarianceEstimator(const EnKGaussianMeasurementCovarianceEstimator &another)
  :GaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void EnKGaussianMeasurementCovarianceEstimator::operator=(const EnKGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  GaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimator(*this));
}

void EnKGaussianMeasurementCovarianceEstimator::make_copy(const EnKGaussianMeasurementCovarianceEstimator &another)
{
}

MeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new EnKGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new EnKGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}
