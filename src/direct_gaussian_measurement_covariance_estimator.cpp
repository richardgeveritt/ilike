#include "direct_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator()
  :GaussianMeasurementCovarianceEstimator()
{
}

DirectGaussianMeasurementCovarianceEstimator::~DirectGaussianMeasurementCovarianceEstimator()
{
}

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another)
  :GaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void DirectGaussianMeasurementCovarianceEstimator::operator=(const DirectGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  GaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimator(*this));
}

void DirectGaussianMeasurementCovarianceEstimator::make_copy(const DirectGaussianMeasurementCovarianceEstimator &another)
{
  //this->measurement_kernel = another.measurement_kernel;
  //this->kernel = another.kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  this->measurement_noise_functions = another.measurement_noise_functions;
  this->transform_function = another.transform_function;
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
}

MeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new DirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new DirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

/*
 arma::mat DirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
 {
 return this->kernel.get_covariance(this->measurement_variables);
 }
 */
