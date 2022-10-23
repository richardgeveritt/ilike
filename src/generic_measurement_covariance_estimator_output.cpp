#include "generic_measurement_covariance_estimator_output.h"
#include "generic_measurement_covariance_estimator.h"

GenericMeasurementCovarianceEstimatorOutput::GenericMeasurementCovarianceEstimatorOutput()
  :MeasurementCovarianceEstimatorOutput()
{
}

GenericMeasurementCovarianceEstimatorOutput::~GenericMeasurementCovarianceEstimatorOutput()
{
}

GenericMeasurementCovarianceEstimatorOutput::GenericMeasurementCovarianceEstimatorOutput(GenericMeasurementCovarianceEstimator* generic_estimator_in)
{
  this->generic_estimator = generic_estimator_in;
}

GenericMeasurementCovarianceEstimatorOutput::GenericMeasurementCovarianceEstimatorOutput(const GenericMeasurementCovarianceEstimatorOutput &another)
  :MeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimatorOutput::operator=(const GenericMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;

  MeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimatorOutput::make_copy(const GenericMeasurementCovarianceEstimatorOutput &another)
{
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  this->generic_estimator = another.generic_estimator;
  //this->simulated_measurement = another.simulated_measurement;
  this->measurement_state = another.measurement_state;
  this->random_shift = another.random_shift;
}

MeasurementCovarianceEstimatorOutput* GenericMeasurementCovarianceEstimatorOutput::duplicate() const
{
  return( new GenericMeasurementCovarianceEstimatorOutput(*this));
}

void GenericMeasurementCovarianceEstimatorOutput::simulate(const Parameters &current_state)
{
  this->measurement_state = this->generic_estimator->simulate(current_state);
  this->random_shift = this->generic_estimator->gaussian_simulate();
}

arma::rowvec GenericMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance() const
{
  return arma::conv_to<arma::rowvec>::from(this->measurement_state);
}

arma::colvec GenericMeasurementCovarianceEstimatorOutput::get_shift(double inverse_incremental_temperature) const
{
  return this->measurement_state + arma::chol((inverse_incremental_temperature-1.0)*this->generic_estimator->Cygivenx)*this->random_shift;
}

arma::colvec GenericMeasurementCovarianceEstimatorOutput::get_deterministic_shift() const
{
  return this->measurement_state;
}

arma::mat GenericMeasurementCovarianceEstimatorOutput::get_kalman_gain(const arma::mat &Cxy,
                                                                       const arma::mat &Cyy,
                                                                       double inverse_incremental_temperature)
{
  return Cxy*(Cyy + (inverse_incremental_temperature-1.0)*this->generic_estimator->get_Cygivenx()).i();
}

arma::mat GenericMeasurementCovarianceEstimatorOutput::get_adjustment(const arma::mat &Zf,
                                                                      const arma::mat &Ginv,
                                                                      const arma::mat &Ftranspose,
                                                                      const arma::mat &V,
                                                                      double inverse_incremental_temperature)
{
  // follows https://arxiv.org/abs/2006.02941
  arma::mat for_eig = V*((inverse_incremental_temperature-1.0)*this->generic_estimator->get_Cygivenx())*V.t();
  
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

double GenericMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                       const Parameters &conditioned_on_parameters)
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                    const Parameters &conditioned_on_parameters)
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::evaluate_likelihood()
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::evaluate_likelihood(const Parameters &conditioned_on_parameters)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood()
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}

double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)
{
  return dmvnorm(*this->estimator->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}
