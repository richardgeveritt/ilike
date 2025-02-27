#include "generic_measurement_covariance_estimator_output.h"
#include "generic_measurement_covariance_estimator.h"
#include "transform.h"

namespace ilike
{
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

void GenericMeasurementCovarianceEstimatorOutput::specific_simulate(const Parameters &current_state)
{
  std::shared_ptr<Transform> summary_statistics = this->generic_estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    this->measurement_state = this->generic_estimator->simulate(current_state).get_colvec(this->generic_estimator->measurement_variables);
  }
  else
  {
    this->measurement_state = summary_statistics->transform(this->generic_estimator->simulate(current_state)).get_colvec(this->generic_estimator->measurement_variables);
  }
  
  this->random_shift = this->generic_estimator->gaussian_simulate();
}

void GenericMeasurementCovarianceEstimatorOutput::subsample_specific_simulate(const Parameters &current_state)
{
  std::shared_ptr<Transform> summary_statistics = this->generic_estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    this->measurement_state = this->generic_estimator->simulate(current_state).get_colvec(this->generic_estimator->measurement_variables);
  }
  else
  {
    this->measurement_state = summary_statistics->transform(this->generic_estimator->simulate(current_state)).get_colvec(this->generic_estimator->measurement_variables);
  }
  this->random_shift = this->generic_estimator->gaussian_simulate();
}

arma::rowvec GenericMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance() const
{
  return arma::conv_to<arma::rowvec>::from(this->measurement_state);
}

arma::colvec GenericMeasurementCovarianceEstimatorOutput::get_shift(double inverse_incremental_temperature) const
{
  if ((inverse_incremental_temperature-1.0)==0.0)
    return this->measurement_state;
  else
    return this->measurement_state + arma::chol((inverse_incremental_temperature-1.0)*this->generic_estimator->Cygivenx)*this->random_shift;
}

arma::colvec GenericMeasurementCovarianceEstimatorOutput::get_deterministic_shift() const
{
  return this->measurement_state;
}

double GenericMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                       const arma::mat &inv_sigma_precomp,
                                                                                       double log_det_precomp) const
{
  double c = (double(inv_sigma_precomp.n_rows)/2.0)*(1.0/inverse_incremental_temperature)*log(inverse_incremental_temperature) + (double(inv_sigma_precomp.n_rows)/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log_det_precomp;
  return c + dmvnorm_using_precomp(*this->generic_estimator->get_measurement_pointer(),
                               this->measurement_state,
                               inv_sigma_precomp,
                               log_det_precomp);
}

/*
 double GenericMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
 }
 */

double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                                 const arma::mat &inv_sigma_precomp,
                                                                                                 double log_det_precomp) const
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  double c = (double(inv_sigma_precomp.n_rows)/2.0)*(1.0/inverse_incremental_temperature)*log(inverse_incremental_temperature) + (double(inv_sigma_precomp.n_rows)/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log_det_precomp;
  return c + dmvnorm(*this->generic_estimator->get_measurement_pointer(),
                 this->measurement_state,
                 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
}

/*
 double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 inverse_incremental_temperature*this->generic_estimator->Cygivenx);
 }
 */

double GenericMeasurementCovarianceEstimatorOutput::evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}

/*
 double GenericMeasurementCovarianceEstimatorOutput::evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 this->generic_estimator->Cygivenx);
 }
 */

double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 this->measurement_state,
                 this->generic_estimator->Cygivenx);
}

/*
 double GenericMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 this->generic_estimator->Cygivenx);
 }
 */

MeasurementCovarianceEstimator* GenericMeasurementCovarianceEstimatorOutput::get_estimator()
{
  return this->generic_estimator;
}

void GenericMeasurementCovarianceEstimatorOutput::write_to_file(const std::string &directory_name,
                                                                const std::string &index)
{
  
}

void GenericMeasurementCovarianceEstimatorOutput::close_ofstreams()
{
  
}
}
