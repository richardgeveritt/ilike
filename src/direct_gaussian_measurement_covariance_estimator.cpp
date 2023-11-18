#include "direct_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator()
  :GaussianMeasurementCovarianceEstimator()
{
}

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                           size_t* seed_in,
                                                                                           Data* data_in,
                                                                                           std::shared_ptr<Transform> transform_in,
                                                                                           std::shared_ptr<Transform> summary_statistics_in)
: GaussianMeasurementCovarianceEstimator(rng_in,
                                         seed_in,
                                         data_in,
                                         transform_in,
                                         summary_statistics_in)
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

/*
MeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* DirectGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new DirectGaussianMeasurementCovarianceEstimator(*this));
}
*/

void DirectGaussianMeasurementCovarianceEstimator::make_copy(const DirectGaussianMeasurementCovarianceEstimator &another)
{
  this->kernel = another.kernel;
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

void DirectGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

/*
void DirectGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  //Data dummy_data = this->transform_function(Parameters());
  //this->measurement_variables = dummy_data.get_vector_variables();
  Data dummy_data = this->transform_function->transform(Parameters());
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(dummy_data[this->measurement_variables[i]].n_elem));
    this->kernel.set_covariance(this->measurement_variables[i],
                                this->measurement_noise_functions[i](Parameters()));
  }
}

void DirectGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data = this->transform_function->transform(conditioned_on_parameters);
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(dummy_data[this->measurement_variables[i]].n_elem));
    this->kernel.set_covariance(this->measurement_variables[i],
                                this->measurement_noise_functions[i](conditioned_on_parameters));
  }
  
  //this->measurement_variables = dummy_data.get_vector_variables();
}
*/

arma::mat DirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::mat DirectGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

/*
arma::mat DirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
return this->kernel.get_covariance(this->measurement_variables);
}
*/

/*
void DirectGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    //this->conditioned_on_parameters = conditioned_on_parameters_in;
    for (size_t i=0;
         i<this->measurement_variables.size();
         ++i)
    {
      this->get_kernel().set_covariance(this->measurement_variables[i],
                                        this->measurement_noise_functions[i](conditioned_on_parameters_in));
    }
  }
}
*/

arma::colvec DirectGaussianMeasurementCovarianceEstimator::get_measurement_state(const Parameters &parameters) const
{
  return this->get_measurement_state_parameters(parameters).get_colvec(this->measurement_variables);
}
