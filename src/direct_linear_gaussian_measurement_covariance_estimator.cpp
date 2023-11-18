#include "direct_linear_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

DirectLinearGaussianMeasurementCovarianceEstimator::DirectLinearGaussianMeasurementCovarianceEstimator()
  :DirectGaussianMeasurementCovarianceEstimator()
{
}

DirectLinearGaussianMeasurementCovarianceEstimator::DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                       size_t* seed_in,
                                                                                                       Data* data_in,
                                                                                                       std::shared_ptr<Transform> transform_in,
                                                                                                       std::shared_ptr<Transform> summary_statistics_in,
                                                                                                       const arma::mat &measurement_matrix_in,
                                                                                                       const arma::mat &measurement_covariance_in,
                                                                                                       const std::string &measurement_variable_in,
                                                                                                       const std::string &state_variable_in)
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
  this->As.push_back(measurement_matrix_in);
  this->state_variable = state_variable_in;
}

DirectLinearGaussianMeasurementCovarianceEstimator::DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                       size_t* seed_in,
                                                                                                       Data* data_in,
                                                                                                       std::shared_ptr<Transform> transform_in,
                                                                                                       std::shared_ptr<Transform> summary_statistics_in,
                                                                                                       GetMatrixPtr measurement_matrix_in,
                                                                                                       GetMatrixPtr measurement_covariance_in,
                                                                                                       const std::string &measurement_variable_in,
                                                                                                       const std::string &state_variable_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  this->set_using_parameters = true;
  this->measurement_variables.push_back(measurement_variable_in);
  this->measurement_noise_functions.push_back(measurement_covariance_in);
  this->A_functions.push_back(measurement_matrix_in);
  this->state_variable = state_variable_in;
}

DirectLinearGaussianMeasurementCovarianceEstimator::~DirectLinearGaussianMeasurementCovarianceEstimator()
{
}

DirectLinearGaussianMeasurementCovarianceEstimator::DirectLinearGaussianMeasurementCovarianceEstimator(const DirectLinearGaussianMeasurementCovarianceEstimator &another)
  :DirectGaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void DirectLinearGaussianMeasurementCovarianceEstimator::operator=(const DirectLinearGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  DirectGaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* DirectLinearGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new DirectLinearGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* DirectLinearGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new DirectLinearGaussianMeasurementCovarianceEstimator(*this));
}

void DirectLinearGaussianMeasurementCovarianceEstimator::make_copy(const DirectLinearGaussianMeasurementCovarianceEstimator &another)
{
  //this->measurement_kernel = another.measurement_kernel;
  //this->kernel = another.kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  this->measurement_noise_functions = another.measurement_noise_functions;
  this->As = another.As;
  this->A_functions = another.A_functions;
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  this->state_variable = another.state_variable;
}

/*
MeasurementCovarianceEstimatorOutput* DirectLinearGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new DirectLinearGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* DirectLinearGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new DirectLinearGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}
*/

/*
void DirectLinearGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectLinearGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}
*/

void DirectLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  //Data dummy_data = this->transform_function(Parameters());
  //this->measurement_variables = dummy_data.get_vector_variables();
  //Data dummy_data = this->transform_function->transform(Parameters());
  
  bool set_A = false;
  
  if (As.size()!=this->measurement_variables.size())
  {
    set_A = true;
    this->As.reserve(this->measurement_variables.size());
  }
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    arma::mat cov = this->measurement_noise_functions[i](Parameters());
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(cov.n_rows));
    this->kernel.set_covariance(this->measurement_variables[i],
                                cov);
    if (set_A)
    {
      this->As.push_back(this->A_functions[i](Parameters()));
    }
  }
}

void DirectLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  bool set_A = false;
  
  if (As.size()!=this->measurement_variables.size())
  {
    set_A = true;
    this->As.reserve(this->measurement_variables.size());
  }
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    arma::mat dummy_data = this->A_functions[i](conditioned_on_parameters)*conditioned_on_parameters[state_variable];
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(dummy_data.n_elem));
    this->kernel.set_covariance(this->measurement_variables[i],
                                this->measurement_noise_functions[i](conditioned_on_parameters));
    if (set_A)
    {
      this->As.push_back(this->A_functions[i](conditioned_on_parameters));
    }
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
arma::mat DirectLinearGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::mat DirectLinearGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}
*/

/*
arma::mat DirectLinearGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
return this->kernel.get_covariance(this->measurement_variables);
}
*/

void DirectLinearGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
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
      this->As[i] = this->A_functions[i](conditioned_on_parameters_in);
    }
  }
}

Parameters DirectLinearGaussianMeasurementCovarianceEstimator::get_measurement_state_parameters(const Parameters &parameters) const
{
  Parameters output;
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    output[measurement_variables[i]] = this->As[i]*parameters[state_variable];
  }
  return output;
}

/*
GaussianIndependentProposalKernel DirectLinearGaussianMeasurementCovarianceEstimator::get_kernel() const
{
  return this->kernel;
}
*/
