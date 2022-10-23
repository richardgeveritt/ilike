#include "gaussian_measurement_covariance_estimator_output.h"
#include "gaussian_measurement_covariance_estimator.h"
#include "distributions.h"

GaussianMeasurementCovarianceEstimatorOutput::GaussianMeasurementCovarianceEstimatorOutput()
  :MeasurementCovarianceEstimatorOutput()
{
}

GaussianMeasurementCovarianceEstimatorOutput::~GaussianMeasurementCovarianceEstimatorOutput()
{
}

GaussianMeasurementCovarianceEstimatorOutput::GaussianMeasurementCovarianceEstimatorOutput(GaussianMeasurementCovarianceEstimator* gaussian_estimator_in)
:MeasurementCovarianceEstimatorOutput(gaussian_estimator_in)
{
}

GaussianMeasurementCovarianceEstimatorOutput::GaussianMeasurementCovarianceEstimatorOutput(const GaussianMeasurementCovarianceEstimatorOutput &another)
  :MeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimatorOutput::operator=(const GaussianMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;

  MeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimatorOutput::make_copy(const GaussianMeasurementCovarianceEstimatorOutput &another)
{
  this->measurement_state = another.measurement_state;
  this->random_shift = another.random_shift;
}

arma::rowvec GaussianMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance() const
{
  return arma::conv_to<arma::rowvec>::from(this->measurement_state);
}

arma::mat GaussianMeasurementCovarianceEstimatorOutput::get_kalman_gain(const arma::mat &Cxy,
                                                                        const arma::mat &Cyy,
                                                                        double inverse_incremental_temperature)
{
  return Cxy*(Cyy + inverse_incremental_temperature*this->get_measurement_covariance()).i();
}

arma::colvec GaussianMeasurementCovarianceEstimatorOutput::get_shift(double inverse_incremental_temperature) const
{
  return this->measurement_state + sqrt(inverse_incremental_temperature)*this->random_shift;
}

arma::colvec GaussianMeasurementCovarianceEstimatorOutput::get_deterministic_shift() const
{
  return this->measurement_state;
}

arma::mat GaussianMeasurementCovarianceEstimatorOutput::get_adjustment(const arma::mat &Zf,
                                                                       const arma::mat &Ginv,
                                                                       const arma::mat &Ftranspose,
                                                                       const arma::mat &V,
                                                                       double inverse_incremental_temperature)
{
  // follows https://arxiv.org/abs/2006.02941
  arma::mat for_eig = V*(inverse_incremental_temperature*this->get_measurement_covariance())*V.t();
  
  arma::mat C;
  arma::vec diagGamma;
  arma::mat Ctrans;
  arma::svd_econ(C,diagGamma,Ctrans,for_eig);
  
  arma::mat Gamma;
  Gamma.diag() = diagGamma;
  arma::mat I;
  I.eye( arma::size(Gamma) );
  
  return Zf*C*arma::sqrtmat_sympd(arma::inv_sympd(I+Gamma))*Ginv*Ftranspose;
}

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                        const Parameters &conditioned_on_parameters)
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                                  const Parameters &conditioned_on_parameters)
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood()
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood(const Parameters &conditioned_on_parameters)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood()
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_measurement_covariance());
}

double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_measurement_covariance());
}
