#include "enk_gaussian_measurement_covariance_estimator.h"
#include "enk_gaussian_measurement_covariance_estimator_output.h"
#include "parameters.h"

namespace ilike
{
EnKGaussianMeasurementCovarianceEstimator::EnKGaussianMeasurementCovarianceEstimator()
:GaussianMeasurementCovarianceEstimator()
{
}
EnKGaussianMeasurementCovarianceEstimator::~EnKGaussianMeasurementCovarianceEstimator()
{
}

EnKGaussianMeasurementCovarianceEstimator::EnKGaussianMeasurementCovarianceEstimator(const EnKGaussianMeasurementCovarianceEstimator &another)
:GaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void EnKGaussianMeasurementCovarianceEstimator::operator=(const EnKGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;
  
  GaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

MeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimator::duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimator(*this));
}

GaussianMeasurementCovarianceEstimator* EnKGaussianMeasurementCovarianceEstimator::gaussian_duplicate() const
{
  return( new EnKGaussianMeasurementCovarianceEstimator(*this));
}

void EnKGaussianMeasurementCovarianceEstimator::make_copy(const EnKGaussianMeasurementCovarianceEstimator &another)
{
}

MeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new EnKGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* EnKGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new EnKGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

void EnKGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void EnKGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

void EnKGaussianMeasurementCovarianceEstimator::setup_measurement_variables()
{
  Data dummy_data;
  this->measurement_variables = dummy_data.get_vector_variables();
}

void EnKGaussianMeasurementCovarianceEstimator::setup_measurement_variables(const Parameters &conditioned_on_parameters)
{
  Data dummy_data;
  this->measurement_variables = dummy_data.get_vector_variables();
}

arma::mat EnKGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return arma::mat();
}

arma::mat EnKGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return arma::mat();
}

arma::mat EnKGaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
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
  arma::eig_sym(diagD,U,for_eig);
  arma::mat Dsqrt(diagD.n_elem,diagD.n_elem);
  Dsqrt.diag() = arma::sqrt(diagD);
  
  return P*Dhathalf*U*Dsqrt*arma::pinv(Dhathalf)*P.t();
  
}

arma::mat EnKGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment(const arma::mat &Cxy,
                                                                         const arma::mat &Cyy,
                                                                         double inverse_incremental_temperature) const
{
  Rcpp::stop("EnKGaussianMeasurementCovarianceEstimator::get_sqrt_adjustment - not yet implemented.");
  
}
}
