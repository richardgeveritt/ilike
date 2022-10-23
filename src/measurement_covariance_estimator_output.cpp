#include "measurement_covariance_estimator_output.h"
#include "measurement_covariance_estimator.h"

MeasurementCovarianceEstimatorOutput::MeasurementCovarianceEstimatorOutput()
{
}

MeasurementCovarianceEstimatorOutput::MeasurementCovarianceEstimatorOutput(MeasurementCovarianceEstimator* estimator_in)
{
  this->estimator = estimator_in;
}

MeasurementCovarianceEstimatorOutput::~MeasurementCovarianceEstimatorOutput()
{
}

MeasurementCovarianceEstimatorOutput::MeasurementCovarianceEstimatorOutput(const MeasurementCovarianceEstimatorOutput &another)
{
  this->make_copy(another);
}

void MeasurementCovarianceEstimatorOutput::operator=(const MeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void MeasurementCovarianceEstimatorOutput::make_copy(const MeasurementCovarianceEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->set_using_parameters = another.set_using_parameters;
}

arma::colvec* MeasurementCovarianceEstimatorOutput::get_measurement() const
{
  return this->estimator->get_measurement_pointer();
}
