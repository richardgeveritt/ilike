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
                                                                                           std::shared_ptr<Transform> summary_statistics_in,
                                                                                           std::shared_ptr<Transform> transform_function_in,
                                                                                           const std::vector<std::string> &measurement_variables_in,
                                                                                           const std::vector<arma::mat> &measurement_noises_in)
: GaussianMeasurementCovarianceEstimator(rng_in,
                                         seed_in,
                                         data_in,
                                         transform_in,
                                         summary_statistics_in)
{
  this->set_using_parameters = false;
  this->measurement_variables = measurement_variables_in;
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(measurement_noises_in[i].n_rows));
    this->kernel.set_covariance(this->measurement_variables[i],
                                measurement_noises_in[i]);
  }
  this->transform_function = transform_function_in;
}

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                           size_t* seed_in,
                                                                                           Data* data_in,
                                                                                           std::shared_ptr<Transform> transform_in,
                                                                                           std::shared_ptr<Transform> summary_statistics_in,
                                                                                           std::shared_ptr<Transform> transform_function_in,
                                                                                           const std::vector<std::string> &measurement_variables_in,
                                                                                           const std::vector<GetMeasurementMatrixPtr> &measurement_noise_functions_in)
: GaussianMeasurementCovarianceEstimator(rng_in,
                                         seed_in,
                                         data_in,
                                         transform_in,
                                         summary_statistics_in)
{
  this->set_using_parameters = true;
  this->measurement_variables = measurement_variables_in;
  this->transform_function = transform_function_in;
  this->measurement_noise_functions = measurement_noise_functions_in;
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
  this->kernel = another.kernel;
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

void DirectGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

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
  
  /*
  for (auto i=dummy_data.vector_begin();
       i!=dummy_data.vector_end();
       ++i)
  {
    this->kernel.set_mean(i->first,arma::colvec(i->second.first->n_elem));
  }
  */
  //std::cout << dummy_data << std::endl;
  //this->measurement_variables = dummy_data.get_vector_variables();
}

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

void DirectGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    //this->conditioned_on_parameters = conditioned_on_parameters_in;
    for (size_t i=0;
         i<this->measurement_variables.size();
         ++i)
    {
      this->kernel.set_covariance(this->measurement_variables[i],
                                  this->measurement_noise_functions[i](conditioned_on_parameters_in));
    }
  }
}
