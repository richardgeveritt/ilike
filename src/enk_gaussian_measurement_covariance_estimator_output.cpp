#include "enk_gaussian_measurement_covariance_estimator_output.h"
#include "enk_gaussian_measurement_covariance_estimator.h"

namespace ilike
{
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

void EnKGaussianMeasurementCovarianceEstimatorOutput::specific_simulate(const Parameters &parameters)
{
  
}

void EnKGaussianMeasurementCovarianceEstimatorOutput::subsample_specific_simulate(const Parameters &parameters)
{
  
}

MeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimatorOutput::get_estimator()
{
  return this->enk_estimator;
}

GaussianMeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimatorOutput::get_gaussian_estimator()
{
  return this->enk_estimator;
}

const GaussianMeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimatorOutput::get_gaussian_estimator() const
{
  return this->enk_estimator;
}
}
