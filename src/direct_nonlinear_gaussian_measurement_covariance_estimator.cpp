#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

namespace ilike
{
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
//const std::string &measurement_variable_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               {measurement_variable_in})
{
  this->set_using_parameters = false;
  
  //if (this->measurement_variables.size()!=1)
  //  Rcpp::stop("DirectGaussianMeasurementCovarianceEstimator - For this EnK likelihood we can only have one variable present in the data.");
  
  this->kernel.set_mean(this->measurement_variables[0],
                        arma::colvec(measurement_covariance_in.n_rows));
  this->kernel.set_covariance(this->measurement_variables[0],
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
//const std::vector<std::string> &measurement_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               measurement_variables_in)
{
  this->set_using_parameters = false;
  
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
//const std::string &measurement_variable_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               {measurement_variable_in})
{
  this->set_using_parameters = true;

  //if (this->measurement_variables.size()!=1)
  //  Rcpp::stop("DirectGaussianMeasurementCovarianceEstimator - For this EnK likelihood we can only have one variable present in the data.");
  
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
//const std::vector<std::string> &measurement_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               measurement_variables_in)
{
  this->set_using_parameters = true;
  
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
  this->measurement_noise_functions = another.measurement_noise_functions;
  this->transform_function = another.transform_function;
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  
  if (this->set_using_parameters)
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
  }
  else
  {
  }
  
}

void DirectNonLinearGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
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

arma::mat DirectNonLinearGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                                const arma::mat &Dhathalf,
                                                                                const arma::mat &P,
                                                                                const arma::mat &Vtranspose,
                                                                                const arma::mat &Yhat,
                                                                                double inverse_incremental_temperature) const
{
  
  // follows https://arxiv.org/abs/2006.02941
  
  arma::mat for_eig = adjustment_inner_matrix(Vtranspose, Yhat, inverse_incremental_temperature*this->get_measurement_covariance());
  
  for_eig = (for_eig+for_eig.t())/2.0;
  
  arma::mat U;
  arma::vec diagD;
  arma::eig_sym(diagD,U,for_eig);
  arma::mat Dsqrt(diagD.n_elem,diagD.n_elem);
  Dsqrt.diag() = arma::sqrt(diagD);
  
  return P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t();
  
}

arma::mat DirectNonLinearGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment(const arma::mat &Cxy,
                                                                                     const arma::mat &Cyy,
                                                                                     double inverse_incremental_temperature) const
{
  Rcpp::stop("DirectNonLinearGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment - not yet implemented.");
}
}
