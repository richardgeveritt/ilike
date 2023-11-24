#include "direct_abc_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

DirectABCGaussianMeasurementCovarianceEstimator::DirectABCGaussianMeasurementCovarianceEstimator()
  :DirectGaussianMeasurementCovarianceEstimator()
{
}

DirectABCGaussianMeasurementCovarianceEstimator::DirectABCGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                 size_t* seed_in,
                                                                                                 Data* data_in,
                                                                                                 std::shared_ptr<Transform> transform_in,
                                                                                                 std::shared_ptr<Transform> summary_statistics_in,
                                                                                                 double min_epsilon_in,
                                                                                                 //const std::string &tempering_variable_in,
                                                                                                 const std::string &scale_variable_in,
                                                                                                 const std::vector<std::string> &measurement_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in)
{
  if (min_epsilon_in<=0.0)
  {
    Rcpp::stop("DirectABCGaussianMeasurementCovarianceEstimator - min_epsilon must be positive.");
  }
  
  this->set_using_parameters = true;
  this->measurement_variables = measurement_variables_in;
  this->min_epsilon = min_epsilon_in;
  //this->tempering_variable = tempering_variable_in;
  this->scale_variable = scale_variable_in;
}

DirectABCGaussianMeasurementCovarianceEstimator::~DirectABCGaussianMeasurementCovarianceEstimator()
{
}

DirectABCGaussianMeasurementCovarianceEstimator::DirectABCGaussianMeasurementCovarianceEstimator(const DirectABCGaussianMeasurementCovarianceEstimator &another)
  :DirectGaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void DirectABCGaussianMeasurementCovarianceEstimator::operator=(const DirectABCGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  DirectGaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* DirectABCGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new DirectABCGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* DirectABCGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new DirectABCGaussianMeasurementCovarianceEstimator(*this));
}

void DirectABCGaussianMeasurementCovarianceEstimator::make_copy(const DirectABCGaussianMeasurementCovarianceEstimator &another)
{
  //this->measurement_kernel = another.measurement_kernel;
  //this->kernel = another.kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  //this->tempering_variable = another.tempering_variable;
  this->scale_variable = another.scale_variable;
  this->min_epsilon = another.min_epsilon;
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
}

/*
MeasurementCovarianceEstimatorOutput* DirectABCGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new DirectABCGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* DirectABCGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new DirectABCGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}
*/

/*
void DirectABCGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectABCGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}
*/

void DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  //Data dummy_data = this->transform_function(Parameters());
  //this->measurement_variables = dummy_data.get_vector_variables();
  
  Rcpp::stop("DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables() - not yet written");
}

void DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data;
  
  //double tempering = conditioned_on_parameters[this->tempering_variable][0];
  
  if (this->summary_statistics==NULL)
  {
    dummy_data = conditioned_on_parameters;
  }
  else
  {
    dummy_data = this->summary_statistics->transform(conditioned_on_parameters);
  }
  
  arma::colvec total_scale;
  if (this->scale_variable!="")
    total_scale = conditioned_on_parameters[this->scale_variable];
  
  size_t current_index = 0;
  
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    size_t current_dummy_data_size = dummy_data[this->measurement_variables[i]].n_elem;
    
    this->kernel.set_mean(this->measurement_variables[i],
                          arma::colvec(dummy_data[this->measurement_variables[i]].n_elem));
    
    arma::colvec scale;
    if (this->scale_variable!="")
    {
      scale = total_scale.rows(current_index,current_index + current_dummy_data_size - 1);
    }
    else
    {
      scale = arma::colvec(current_dummy_data_size);
      scale.fill(1.0);
    }
    
    current_index = current_index + current_dummy_data_size;

    arma::mat cov(scale.n_elem,scale.n_elem,arma::fill::zeros);
    //arma::vec for_diag = (1.0/tempering)*arma::pow(this->min_epsilon*scale,2.0);
    arma::vec for_diag = arma::pow(this->min_epsilon*scale,2.0);
    cov.diag() = for_diag;
    
    this->kernel.set_covariance(this->measurement_variables[i],
                                cov);
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
arma::mat DirectABCGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::mat DirectABCGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}
*/

/*
arma::mat DirectABCGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
return this->kernel.get_covariance(this->measurement_variables);
}
*/

void DirectABCGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    //double tempering = conditioned_on_parameters_in[this->tempering_variable][0];
    
    arma::colvec total_scale;
    if (this->scale_variable!="")
      total_scale = conditioned_on_parameters_in[this->scale_variable];
    
    size_t current_index = 0;
    
    //this->conditioned_on_parameters = conditioned_on_parameters_in;
    for (size_t i=0;
         i<this->measurement_variables.size();
         ++i)
    {
      size_t current_dummy_data_size = this->kernel.get_covariance(measurement_variables[i]).n_rows;
      
      arma::colvec scale;
      if (this->scale_variable!="")
      {
        scale = total_scale.rows(current_index,current_index + current_dummy_data_size - 1);
      }
      else
      {
        scale = arma::colvec(current_dummy_data_size);
        scale.fill(1.0);
      }
      
      current_index = current_index + current_dummy_data_size;
      
      arma::mat cov(scale.n_elem,scale.n_elem,arma::fill::zeros);
      //arma::vec for_diag = (1.0/tempering)*arma::pow(this->min_epsilon*scale,2.0);
      arma::vec for_diag = arma::pow(this->min_epsilon*scale,2.0);
      cov.diag() = for_diag;
      
      this->kernel.set_covariance(this->measurement_variables[i],
                                  cov);

    }
  }
}

Parameters DirectABCGaussianMeasurementCovarianceEstimator::get_measurement_state_parameters(const Parameters &parameters) const
{
  return parameters;
}

/*
GaussianIndependentProposalKernel DirectABCGaussianMeasurementCovarianceEstimator::get_kernel() const
{
  return this->kernel;
}
*/
