#include "direct_linear_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

namespace ilike
{
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
                                                                                                       const std::vector<std::string> &state_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               {measurement_variable_in})
{
  this->set_using_parameters = false;
  //this->measurement_variables.push_back(measurement_variable_in); // reinstate if we don't automatically get from data
  
  //if (this->measurement_variables.size()!=1)
  //  Rcpp::stop("DirectGaussianMeasurementCovarianceEstimator - for this EnK likelihood we can only .");
  
  this->kernel.set_mean(this->measurement_variables[0],
                        arma::colvec(measurement_covariance_in.n_rows));
  this->kernel.set_covariance(this->measurement_variables[0],
                              measurement_covariance_in);
  this->As.push_back(measurement_matrix_in);
  this->state_variables = state_variables_in;
}

DirectLinearGaussianMeasurementCovarianceEstimator::DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                                       size_t* seed_in,
                                                                                                       Data* data_in,
                                                                                                       std::shared_ptr<Transform> transform_in,
                                                                                                       std::shared_ptr<Transform> summary_statistics_in,
                                                                                                       GetMatrixPtr measurement_matrix_in,
                                                                                                       GetMatrixPtr measurement_covariance_in,
                                                                                                       const std::string &measurement_variable_in,
                                                                                                       const std::vector<std::string> &state_variables_in)
: DirectGaussianMeasurementCovarianceEstimator(rng_in,
                                               seed_in,
                                               data_in,
                                               transform_in,
                                               summary_statistics_in,
                                               {measurement_variable_in})
{
  this->set_using_parameters = true;
  
  if (this->measurement_variables.size()!=1)
    Rcpp::stop("DirectGaussianMeasurementCovarianceEstimator - for this EnK likelihood we can only .");
  
  this->measurement_noise_functions.push_back(measurement_covariance_in);
  this->A_functions.push_back(measurement_matrix_in);
  this->state_variables = state_variables_in;
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
  this->measurement_noise_functions = another.measurement_noise_functions;
  this->As = another.As;
  this->A_functions = another.A_functions;
  this->state_variables = another.state_variables;
}

void DirectLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  
}

void DirectLinearGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  if (this->set_using_parameters==true)
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
      arma::mat dummy_data = this->A_functions[i](conditioned_on_parameters)*conditioned_on_parameters[this->state_variables];
      this->kernel.set_mean(this->measurement_variables[i],
                            arma::colvec(dummy_data.n_elem));
      this->kernel.set_covariance(this->measurement_variables[i],
                                  this->measurement_noise_functions[i](conditioned_on_parameters));
      if (set_A)
      {
        this->As.push_back(this->A_functions[i](conditioned_on_parameters));
      }
    }
  }
  else
  {
  }
  
}

void DirectLinearGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
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
    output[measurement_variables[i]] = this->As[i]*parameters[this->state_variables];
  }
  return output;
}

// (C264) the stacked measurement state, as
// get_measurement_state_parameters(parameters).get_colvec(measurement_variables)
// gives it (each A_i x vectorised, in order), without building a Parameters
// for every ensemble member at every step
arma::colvec DirectLinearGaussianMeasurementCovarianceEstimator::get_measurement_state(const Parameters &parameters) const
{
  if (this->measurement_variables.size() == 1)
    return arma::vectorise(arma::mat(this->As[0]*parameters[this->state_variables]));
  arma::colvec stacked;
  for (size_t i=0;
       i<this->measurement_variables.size();
       ++i)
  {
    arma::colvec part = arma::vectorise(arma::mat(this->As[i]*parameters[this->state_variables]));
    if (i == 0)
      stacked = part;
    else
      stacked = arma::join_cols(stacked, part);
  }
  return stacked;
}

arma::mat DirectLinearGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
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

arma::mat DirectLinearGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment(const arma::mat &Cxy,
                                                                                  const arma::mat &Cyy,
                                                                                  double inverse_incremental_temperature) const
{
  
  // Symmetric square roots: the square-root update (I - K H) with
  // K = Cxy S^{-1/2} (S^{1/2} + V^{1/2})^{-1} gives the Kalman posterior covariance
  // only for symmetric roots (a Cholesky factor is triangular; B19).
  arma::mat sqrtV = arma::sqrtmat_sympd(inverse_incremental_temperature*this->get_measurement_covariance());
  arma::mat sqrtS = arma::sqrtmat_sympd(Cyy + inverse_incremental_temperature*this->get_measurement_covariance());
  
  arma::mat stacked_H = this->As[0];
  if (this->As.size()>1)
  {
    for (size_t i=1; i<this->As.size(); ++i)
    {
      stacked_H = arma::join_cols(stacked_H, As[i]);
    }
  }
  
  arma::mat K = Cxy * arma::inv_sympd(sqrtS) * arma::inv_sympd(sqrtS + sqrtV);
  
  arma::mat I;
  I.eye(Cxy.n_rows,Cxy.n_rows); // (C278: I - K H is state x state; n_cols, the measurement dimension, worked only when the two were equal)
  
  return I-K*stacked_H;
}
}
