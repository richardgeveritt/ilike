#include "mixed_generic_direct_gaussian_measurement_covariance_estimator_output.h"
#include "mixed_generic_direct_gaussian_measurement_covariance_estimator.h"
#include "transform.h"

namespace ilike
{
MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput()
:MeasurementCovarianceEstimatorOutput()
{
}

MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::~MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput()
{
}

MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(MixedGenericDirectGaussianMeasurementCovarianceEstimator* estimator_in)
{
  this->estimator = estimator_in;
}

MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another)
:MeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;
  
  MeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another)
{
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  this->estimator = another.estimator;
  //this->simulated_measurement = another.simulated_measurement;
  this->likelihood_measurement_state = another.likelihood_measurement_state;
  this->prior_measurement_state = another.prior_measurement_state;
  this->likelihood_random_shift = another.likelihood_random_shift;
  this->prior_random_shift = another.prior_random_shift;
}

MeasurementCovarianceEstimatorOutput* MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::duplicate() const
{
  return( new MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(*this));
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::specific_simulate(const Parameters &parameters)
{
  // do transform on params (could actually leave until later since deterministic, but choose to do now and store result
  std::shared_ptr<Transform> summary_statistics = this->estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    this->likelihood_measurement_state = this->estimator->simulate(parameters).get_colvec(this->estimator->measurement_variables);
  }
  else
  {
    this->likelihood_measurement_state = summary_statistics->transform(this->estimator->simulate(parameters)).get_colvec(this->estimator->measurement_variables);
  }
  
  if (this->estimator->transform_function!=NULL)
  {
    this->prior_measurement_state = this->estimator->transform_function(parameters).get_colvec(this->estimator->prior_measurement_variables);
  }
  else
  {
    this->prior_measurement_state = parameters.get_colvec(this->estimator->prior_measurement_variables);
  }
  
  this->likelihood_random_shift = this->estimator->gaussian_simulate();
  this->estimator->set_parameters(parameters);
  this->prior_random_shift = this->estimator->kernel.independent_simulate(*this->estimator->rng).get_colvec(this->estimator->prior_measurement_variables);
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::subsample_specific_simulate(const Parameters &parameters)
{
  // do transform on params (could actually leave until later since deterministic, but choose to do now and store result
  std::shared_ptr<Transform> summary_statistics = this->estimator->summary_statistics;
  if (summary_statistics==NULL)
  {
    this->likelihood_measurement_state = this->estimator->simulate(parameters).get_colvec(this->estimator->measurement_variables);
  }
  else
  {
    this->likelihood_measurement_state = summary_statistics->transform(this->estimator->simulate(parameters)).get_colvec(this->estimator->measurement_variables);
  }
  
  if (this->estimator->transform_function!=NULL)
  {
    this->prior_measurement_state = this->estimator->transform_function(parameters).get_colvec(this->estimator->prior_measurement_variables);
  }
  else
  {
    this->prior_measurement_state = parameters.get_colvec(this->estimator->prior_measurement_variables);
  }
  
  this->likelihood_random_shift = this->estimator->gaussian_simulate();
  this->estimator->set_parameters(parameters);
  this->prior_random_shift = this->estimator->kernel.independent_simulate(*this->estimator->rng).get_colvec(this->estimator->prior_measurement_variables);
}

arma::rowvec MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance() const
{
  return arma::conv_to<arma::rowvec>::from(arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state));
}

arma::colvec MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::get_shift(double inverse_incremental_temperature) const
{
  arma::colvec blank_measurement = arma::colvec(this->likelihood_measurement_state.n_elem);
  arma::colvec blank_prior = arma::colvec(this->prior_measurement_state.n_elem);
  if ((inverse_incremental_temperature-1.0)==0.0)
  {
    return arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state) + arma::join_cols(blank_measurement,this->prior_random_shift);
  }
  else
  {
    return arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state) + arma::join_cols(arma::chol((inverse_incremental_temperature-1.0)*this->estimator->likelihood_Cygivenx)*this->likelihood_random_shift,blank_prior) + arma::join_cols(blank_measurement,sqrt(inverse_incremental_temperature)*this->prior_random_shift);
  }
}

arma::colvec MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::get_deterministic_shift() const
{
  return arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state);
}

double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                                          const arma::mat &inv_sigma_precomp,
                                                                                                          const double log_det_precomp) const

{
  double c = (double(inv_sigma_precomp.n_rows)/2.0)*(1.0/inverse_incremental_temperature)*log(inverse_incremental_temperature) + (double(inv_sigma_precomp.n_rows)/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log_det_precomp;
  return c + dmvnorm_using_precomp(*this->estimator->get_measurement_pointer(),
                               arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
                               inv_sigma_precomp,
                               log_det_precomp);
}

/*
 double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
 this->estimator->get_measurement_covariance_for_likelihood_ratio(inverse_incremental_temperature));
 }
 */

double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                                                    const arma::mat &inv_sigma_precomp,
                                                                                                                    double log_det_precomp) const
{
  double c = (double(inv_sigma_precomp.n_rows)/2.0)*(1.0/inverse_incremental_temperature)*log(inverse_incremental_temperature) + (double(inv_sigma_precomp.n_rows)/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log_det_precomp;
  // parameters of covariance should already be set at this point, so second argument does nothing
  return c + dmvnorm(*this->estimator->get_measurement_pointer(),
                 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
                 this->estimator->get_measurement_covariance_for_likelihood_ratio(inverse_incremental_temperature));
}

/*
 double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
 this->estimator->get_measurement_covariance_for_likelihood_ratio(inverse_incremental_temperature));
 }
 */

double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
                 this->estimator->get_measurement_covariance_for_likelihood_ratio(1.0));
}

/*
 double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
 this->estimator->get_measurement_covariance_for_likelihood_ratio(1.0));
 }
 */

double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
                 this->estimator->get_measurement_covariance_for_likelihood_ratio(1.0));
}

/*
 double MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 arma::join_cols(this->likelihood_measurement_state,this->prior_measurement_state),
 this->estimator->get_measurement_covariance_for_likelihood_ratio(1.0));
 }
 */

MeasurementCovarianceEstimator* MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::get_estimator()
{
  return this->estimator;
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::write_to_file(const std::string &directory_name,
                                                                                   const std::string &index)
{
  
}

void MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput::close_ofstreams()
{
  
}
}
