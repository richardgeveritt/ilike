#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator()
  :DirectGaussianMeasurementCovarianceEstimator()
{
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                             size_t* seed_in,
                                                                                                             Data* data_in,
                                                                                                             std::shared_ptr<Transform> transform_in,
                                                                                                             std::shared_ptr<Transform> summary_statistics_in,
                                                                                                             std::shared_ptr<Transform> transform_function_in,
                                                                                                             const arma::mat &measurement_covariance_in,
                                                                                                             const std::string &measurement_variable_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  this->set_using_parameters = false;
  this->measurement_variables.push_back(measurement_variable_in);
  this->kernel.set_mean(measurement_variable_in,
                        arma::colvec(measurement_covariance_in.n_rows));
  this->kernel.set_covariance(measurement_variable_in,
                              measurement_covariance_in);
  this->transform_function = transform_function_in;
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                             size_t* seed_in,
                                                                                                             Data* data_in,
                                                                                                             std::shared_ptr<Transform> transform_in,
                                                                                                             std::shared_ptr<Transform> summary_statistics_in,
                                                                                                             std::shared_ptr<Transform> transform_function_in,
                                                                                                             const std::vector<arma::mat> &measurement_covariances_in,
                                                                                                             const std::vector<std::string> &measurement_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  this->set_using_parameters = false;
  
  this->measurement_variables = measurement_variables_in;
  for (size_t i=0; i<this->measurement_variables.size(); ++i)
  {
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(measurement_covariances_in[i].n_rows));
    this->kernel.set_covariance(this->measurement_variables[i],
                                measurement_covariances_in[i]);
  }
  this->transform_function = transform_function_in;
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                             size_t* seed_in,
                                                                                                             Data* data_in,
                                                                                                             std::shared_ptr<Transform> transform_in,
                                                                                                             std::shared_ptr<Transform> summary_statistics_in,
                                                                                                             std::shared_ptr<Transform> transform_function_in,
                                                                                                             GetMatrixPtr measurement_covariance_function_in,
                                                                                                             const std::string &measurement_variable_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  this->set_using_parameters = true;
  this->measurement_variables.push_back(measurement_variable_in);
  this->measurement_noise_functions.push_back(measurement_covariance_function_in);
  this->transform_function = transform_function_in;
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                             size_t* seed_in,
                                                                                                             Data* data_in,
                                                                                                             std::shared_ptr<Transform> transform_in,
                                                                                                             std::shared_ptr<Transform> summary_statistics_in,
                                                                                                             std::shared_ptr<Transform> transform_function_in,
                                                                                                             const std::vector<GetMatrixPtr> &measurement_noise_functions_in,
                                                                                                             const std::vector<std::string> &measurement_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  this->set_using_parameters = true;
  
  this->measurement_variables = measurement_variables_in;
  this->measurement_noise_functions = measurement_noise_functions_in;
  this->transform_function = transform_function_in;
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::~DirectNonLinearGaussianMeasurementCovarianceEstimator()
{
}

DirectNonLinearGaussianMeasurementCovarianceEstimator::DirectNonLinearGaussianMeasurementCovarianceEstimator(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another)
  :DirectGaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::operator=(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  DirectGaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* DirectNonLinearGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new DirectNonLinearGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* DirectNonLinearGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new DirectNonLinearGaussianMeasurementCovarianceEstimator(*this));
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::make_copy(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another)
{
  //this->measurement_kernel = another.measurement_kernel;
  //this->kernel = another.kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  this->measurement_noise_functions = another.measurement_noise_functions;
  this->transform_function = another.transform_function;
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
}

/*
MeasurementCovarianceEstimatorOutput* DirectNonLinearGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new DirectNonLinearGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* DirectNonLinearGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new DirectNonLinearGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}
*/

/*
void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}
*/

void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  //Data dummy_data = this->transform_function(Parameters());
  //this->measurement_variables = dummy_data.get_vector_variables();
  
  Data dummy_data;
  
  if (this->summary_statistics==NULL)
  {
    dummy_data = this->transform_function->transform(Parameters());
  }
  else
  {
    dummy_data = this->summary_statistics->transform(this->transform_function->transform(Parameters()));
  }
  
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

void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data;
  
  if (this->summary_statistics==NULL)
  {
    dummy_data = this->transform_function->transform(conditioned_on_parameters);
  }
  else
  {
    dummy_data = this->summary_statistics->transform(this->transform_function->transform(conditioned_on_parameters));
  }
  
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
  //this->measurement_variables = dummy_data.get_vector_variables();
}

/*
arma::mat DirectNonLinearGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::mat DirectNonLinearGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}
*/

/*
arma::mat DirectNonLinearGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
return this->kernel.get_covariance(this->measurement_variables);
}
*/

void DirectNonLinearGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
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

Parameters DirectNonLinearGaussianMeasurementCovarianceEstimator::get_measurement_state_parameters(const Parameters &parameters) const
{
  if (this->transform_function!=NULL)
  {
    return this->transform_function->transform(parameters);
  }
  else
  {
    return parameters;
  }
}

/*
GaussianIndependentProposalKernel DirectNonLinearGaussianMeasurementCovarianceEstimator::get_kernel() const
{
  return this->kernel;
}
*/
