#include "measurement_covariance_estimator_output.h"
#include "measurement_covariance_estimator.h"
#include "transform.h"

MeasurementCovarianceEstimatorOutput::MeasurementCovarianceEstimatorOutput()
{
  this->write_to_file_flag = true;
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
  this->write_to_file_flag = another.write_to_file_flag;
  //this->estimator = another.estimator;
  //this->set_using_parameters = another.set_using_parameters;
}

void MeasurementCovarianceEstimatorOutput::simulate(const Parameters &parameters)
{
  std::shared_ptr<Transform> transform = this->get_estimator()->transform;
  if (transform==NULL)
    this->specific_simulate(parameters);
  else
    this->specific_simulate(transform->inverse_transform(parameters));
}

void MeasurementCovarianceEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
  std::shared_ptr<Transform> transform = this->get_estimator()->transform;
  if (transform==NULL)
    this->subsample_specific_simulate(parameters);
  else
    this->subsample_specific_simulate(transform->inverse_transform(parameters));
}

arma::colvec* MeasurementCovarianceEstimatorOutput::get_measurement()
{
  return this->get_estimator()->get_measurement_pointer();
}

void MeasurementCovarianceEstimatorOutput::write(const std::string &directory_name)
{
  if (this->write_to_file_flag==true)
    this->write_to_file(directory_name);
}
