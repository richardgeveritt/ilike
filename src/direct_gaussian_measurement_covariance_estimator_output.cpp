#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "direct_gaussian_measurement_covariance_estimator.h"

DirectGaussianMeasurementCovarianceEstimatorOutput::DirectGaussianMeasurementCovarianceEstimatorOutput()
  :GaussianMeasurementCovarianceEstimatorOutput()
{
}

DirectGaussianMeasurementCovarianceEstimatorOutput::DirectGaussianMeasurementCovarianceEstimatorOutput(DirectGaussianMeasurementCovarianceEstimator* direct_estimator_in)
:GaussianMeasurementCovarianceEstimatorOutput(direct_estimator_in)
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
  this->kernel = another.kernel;
  this->direct_estimator = another.direct_estimator;
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::simulate(const Parameters &parameters)
{
  // do transform on params (could actually leave until later since deterministic, but choose to do now and store result
  this->measurement_state = this->direct_estimator->transform_function(parameters).get_vector(this->direct_estimator->measurement_variables);
  
  this->set_parameters(parameters);
  this->random_shift = this->kernel.independent_simulate(*this->direct_estimator->rng).get_vector(this->direct_estimator->measurement_variables);
}

/*
arma::rowvec DirectGaussianMeasurementCovarianceEstimatorOutput::get_measurement_random_shift()
{
  
}
*/

arma::mat DirectGaussianMeasurementCovarianceEstimatorOutput::get_measurement_covariance()
{
  return this->kernel.get_covariance(this->direct_estimator->measurement_variables);
}

void DirectGaussianMeasurementCovarianceEstimatorOutput::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    //this->conditioned_on_parameters = conditioned_on_parameters_in;
    for (size_t i=0;
         i<this->direct_estimator->measurement_variables.size();
         ++i)
    {
      this->kernel.set_covariance(this->direct_estimator->measurement_variables[i],
                                  this->direct_estimator->measurement_noise_functions[i](conditioned_on_parameters_in));
    }
  }
}
