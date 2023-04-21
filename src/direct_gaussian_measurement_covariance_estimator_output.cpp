#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "transform.h"

DirectGaussianMeasurementCovarianceEstimatorOutput::DirectGaussianMeasurementCovarianceEstimatorOutput()
  :GaussianMeasurementCovarianceEstimatorOutput()
{
}

DirectGaussianMeasurementCovarianceEstimatorOutput::DirectGaussianMeasurementCovarianceEstimatorOutput(DirectGaussianMeasurementCovarianceEstimator* direct_estimator_in)
:GaussianMeasurementCovarianceEstimatorOutput()
{
  this->direct_estimator = direct_estimator_in;
}

DirectGaussianMeasurementCovarianceEstimatorOutput::~DirectGaussianMeasurementCovarianceEstimatorOutput()
{
}

DirectGaussianMeasurementCovarianceEstimatorOutput::DirectGaussianMeasurementCovarianceEstimatorOutput(const DirectGaussianMeasurementCovarianceEstimatorOutput &another)
  :GaussianMeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::operator=(const DirectGaussianMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;

  GaussianMeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimatorOutput::duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimatorOutput(*this));
}

GaussianMeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimatorOutput::gaussian_duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimatorOutput(*this));
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::make_copy(const DirectGaussianMeasurementCovarianceEstimatorOutput &another)
{
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  //this->kernel = another.kernel;
  this->direct_estimator = another.direct_estimator;
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::specific_simulate(const Parameters &parameters)
{
  // do transform on params (could actually leave until later since deterministic, but choose to do now and store result
  std::shared_ptr<Transform> summary_statistics = this->direct_estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    if (this->direct_estimator->transform_function!=NULL)
    {
      this->measurement_state = this->direct_estimator->transform_function->transform(parameters).get_colvec(this->direct_estimator->measurement_variables);
    }
    else
    {
      this->measurement_state = parameters.get_colvec(this->direct_estimator->measurement_variables);
    }
  }
  else
  {
    if (this->direct_estimator->transform_function!=NULL)
    {
      this->measurement_state = summary_statistics->transform(this->direct_estimator->transform_function->transform(parameters)).get_colvec(this->direct_estimator->measurement_variables);
    }
    else
    {
      this->measurement_state = summary_statistics->transform(parameters).get_colvec(this->direct_estimator->measurement_variables);
    }
    
  }
  
  this->direct_estimator->set_parameters(parameters);
  this->random_shift = this->direct_estimator->kernel.independent_simulate(*this->direct_estimator->rng).get_colvec(this->direct_estimator->measurement_variables);
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::subsample_specific_simulate(const Parameters &parameters)
{
  // do transform on params (could actually leave until later since deterministic, but choose to do now and store result
  std::shared_ptr<Transform> summary_statistics = this->direct_estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    this->measurement_state = this->direct_estimator->transform_function->transform(parameters).get_colvec(this->direct_estimator->measurement_variables);
  }
  else
  {
    this->measurement_state = summary_statistics->transform(this->direct_estimator->transform_function->transform(parameters)).get_colvec(this->direct_estimator->measurement_variables);
  }
  
  this->direct_estimator->set_parameters(parameters);
  this->random_shift = this->direct_estimator->kernel.independent_simulate(*this->direct_estimator->rng).get_colvec(this->direct_estimator->measurement_variables);
}

/*
arma::rowvec DirectGaussianMeasurementCovarianceEstimatorOutput::get_measurement_random_shift()
{
  
}
*/

MeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimatorOutput::get_estimator()
{
  return this->direct_estimator;
}

GaussianMeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimatorOutput::get_gaussian_estimator()
{
  return this->direct_estimator;
}
