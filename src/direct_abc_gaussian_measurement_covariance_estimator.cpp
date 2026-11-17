#include "direct_abc_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

namespace ilike
{
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
                                               summary_statistics_in,
                                               measurement_variables_in)
{
  if (min_epsilon_in<=0.0)
  {
    Rcpp::stop("DirectABCGaussianMeasurementCovarianceEstimator - min_epsilon must be positive.");
  }
  
  this->set_using_parameters = true;
  this->min_epsilon = min_epsilon_in;
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
  this->scale_variable = another.scale_variable;
  this->min_epsilon = another.min_epsilon;
}

void DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  
  Rcpp::stop("DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables() - not yet written");
}

void DirectABCGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data;

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
    arma::vec for_diag = arma::pow(this->min_epsilon*scale,2.0);
    cov.diag() = for_diag;
    
    this->kernel.set_covariance(this->measurement_variables[i],
                                cov);
  }
  
}

void DirectABCGaussianMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    
    arma::colvec total_scale;
    if (this->scale_variable!="")
      total_scale = conditioned_on_parameters_in[this->scale_variable];
    
    size_t current_index = 0;
    
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

arma::mat DirectABCGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                          const arma::mat &Dhathalf,
                                                                          const arma::mat &P,
                                                                          const arma::mat &Vtranspose,
                                                                          const arma::mat &Yhat,
                                                                          double inverse_incremental_temperature) const
{
  
  // follows https://arxiv.org/abs/2006.02941
  
  arma::mat for_eig = adjustment_inner_matrix(Vtranspose, Yhat, inverse_incremental_temperature*this->get_measurement_covariance());
  
  arma::mat U;
  arma::vec diagD;
  
  for_eig = (for_eig+for_eig.t())/2.0;
  
  arma::eig_sym(diagD,U,for_eig);
  arma::mat Dsqrt(diagD.n_elem,diagD.n_elem);
  Dsqrt.diag() = arma::sqrt(diagD);
  
  return P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t();
  
}

arma::mat DirectABCGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment(const arma::mat &Cxy,
                                                                               const arma::mat &Cyy,
                                                                               double inverse_incremental_temperature) const
{
  // Symmetric square roots: the square-root update (I - K H) with
  // K = Cxy S^{-1/2} (S^{1/2} + V^{1/2})^{-1} gives the Kalman posterior covariance
  // only for symmetric roots (a Cholesky factor is triangular; B19).
  arma::mat sqrtV = arma::sqrtmat_sympd(inverse_incremental_temperature*this->get_measurement_covariance());
  
  arma::colvec abs_eigenvalues = arma::abs(arma::eig_gen(Cyy));
  double biggest_eig = arma::abs(arma::eig_gen(Cyy)).max();
  
  double min_nonzero_eig = 1e50;
  double second_biggest_eig = arma::abs(arma::eig_gen(Cyy)).min();
  for (size_t i=0; i<abs_eigenvalues.n_elem; ++i)
  {
    if ( (abs_eigenvalues[i]<min_nonzero_eig) && (abs_eigenvalues[i]>1e-10) )
    {
      min_nonzero_eig = abs_eigenvalues[i];
    }
    
    if ( (abs_eigenvalues[i]>second_biggest_eig) && (abs_eigenvalues[i]<biggest_eig) )
    {
      second_biggest_eig = abs_eigenvalues[i];
    }
  }
  
  double ratio = arma::abs(arma::eig_gen(Cyy)).max()/min_nonzero_eig;
  
  double scale = ratio/1e14;
  
  arma::mat A;
  A.eye(Cyy.n_rows,Cyy.n_cols);
  A = scale*A;

  arma::mat sqrtS = arma::sqrtmat_sympd(A+Cyy + inverse_incremental_temperature*this->get_measurement_covariance());

  arma::mat stacked_H(Cyy.n_rows,Cxy.n_rows);
  stacked_H.eye();
  
  arma::mat K = Cxy * arma::inv_sympd(sqrtS) * arma::inv_sympd(sqrtS + sqrtV);
  
  arma::mat I;
  I.eye(Cxy.n_rows,Cxy.n_rows);
  
  return I-K*stacked_H;
}
}
