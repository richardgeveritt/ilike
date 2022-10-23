#include "ensemble_kalman_output.h"
#include "ensemble_kalman.h"

EnsembleKalmanOutput::EnsembleKalmanOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  this->estimator = NULL;
}

EnsembleKalmanOutput::EnsembleKalmanOutput(EnsembleKalman* estimator_in,
                                           size_t lag_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  this->lag = lag_in;
}

EnsembleKalmanOutput::~EnsembleKalmanOutput()
{

}

//Copy constructor for the EnsembleKalmanOutput class.
EnsembleKalmanOutput::EnsembleKalmanOutput(const EnsembleKalmanOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void EnsembleKalmanOutput::operator=(const EnsembleKalmanOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* EnsembleKalmanOutput::duplicate() const
{
  return( new EnsembleKalmanOutput(*this));
}

void EnsembleKalmanOutput::make_copy(const EnsembleKalmanOutput &another)
{
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
  this->all_ensembles = another.all_ensembles;
  this->lag = another.lag;
}

Ensemble* EnsembleKalmanOutput::add_ensemble()
{
  size_t num_to_pop_back = std::max<int>(0,this->all_ensembles.size()-lag-1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->all_ensembles.pop_back();
  }
  this->all_ensembles.push_back(Ensemble());
  this->all_ensembles.back().reserve(this->estimator->number_of_ensemble_members);
  return &this->all_ensembles.back();
}

Ensemble EnsembleKalmanOutput::back() const
{
  return this->all_ensembles.back();
}

Ensemble& EnsembleKalmanOutput::back()
{
  return this->all_ensembles.back();
}

void EnsembleKalmanOutput::simulate()
{
  this->estimator->ensemble_kalman_simulate(this);
}

void EnsembleKalmanOutput::simulate(const Parameters &parameters)
{
  this->estimator->ensemble_kalman_simulate(this, parameters);
}

void EnsembleKalmanOutput::evaluate_smcfixed_part(const Parameters &conditioned_on_parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->ensemble_kalman_evaluate(this, conditioned_on_parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->ensemble_kalman_evaluate_smcfixed_part(this, conditioned_on_parameters);
  }
}

void EnsembleKalmanOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(this,conditioned_on_parameters);
  }
  
}

void EnsembleKalmanOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator->ensemble_kalman_simulate(this, parameters);
}

void EnsembleKalmanOutput::subsample_evaluate_smcfixed_part(const Parameters &conditioned_on_parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->ensemble_kalman_subsample_evaluate(this, conditioned_on_parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->ensemble_kalman_subsample_evaluate_smcfixed_part(this, conditioned_on_parameters);
  }
}

void EnsembleKalmanOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(this,conditioned_on_parameters);
  }
  
}

LikelihoodEstimator* EnsembleKalmanOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat EnsembleKalmanOutput::get_gradient_of_log(const std::string &variable,
                                                    const Parameters &x)
{
  throw std::runtime_error("EnsembleKalmanOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat EnsembleKalmanOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                              const Parameters &x)
{
  throw std::runtime_error("EnsembleKalmanOutput::subsample_get_gradient_of_log - not yet implemented.");
}

size_t EnsembleKalmanOutput::number_of_ensemble_kalman_iterations() const
{
  return this->all_ensembles.size();
}


/*
void EnsembleKalmanOutput::set_current_predicted_statistics(const arma::colvec &latest_mean,
                                                          const arma::mat &latest_covariance)
{
  // unsure
  // copied from KF
  this->current_predicted_mean = latest_mean;
  this->current_predicted_covariance = latest_covariance;
}

void EnsembleKalmanOutput::set_current_posterior_statistics(const arma::colvec &latest_mean,
                                                          const arma::mat &latest_covariance)
{
  // unsure
  // copied from KF
  this->current_posterior_mean = latest_mean;
  this->current_posterior_covariance = latest_covariance;
}

void EnsembleKalmanOutput::add_predicted_statistics()
{
  // unsure
  // copied from KF
  size_t num_to_pop_back = std::max<int>(0,this->predicted_means.size()-lag-1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->predicted_means.pop_back();
  }
  this->predicted_means.push_back(this->current_predicted_mean);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->predicted_covariances.pop_back();
  }
  this->predicted_covariances.push_back(this->current_predicted_covariance);
}

void EnsembleKalmanOutput::add_posterior_statistics()
{
  // unsure
  // copied from KF
  size_t num_to_pop_back = std::max<int>(0,this->posterior_means.size()-lag-1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->posterior_means.pop_back();
  }
  this->posterior_means.push_back(this->current_posterior_mean);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->posterior_covariances.pop_back();
  }
  this->posterior_covariances.push_back(this->current_posterior_covariance);
}

void EnsembleKalmanOutput::set_current_predicted_to_be_current_posterior()
{
  // unsure
  // copied from KF
  this->current_predicted_mean = this->current_posterior_mean;
  this->current_predicted_covariance = this->current_posterior_covariance;
}

arma::colvec EnsembleKalmanOutput::predicted_mean_back() const
{
  return this->predicted_means.back();
}

arma::colvec EnsembleKalmanOutput::posterior_mean_back() const
{
  return this->posterior_means.back();
}

arma::mat EnsembleKalmanOutput::predicted_covariance_back() const
{
  return this->predicted_covariances.back();
}

arma::mat EnsembleKalmanOutput::posterior_covariance_back() const
{
  return this->posterior_covariances.back();
}

size_t EnsembleKalmanOutput::predicted_size() const
{
  return this->predicted_means.size();
}

void EnsembleKalmanOutput::print(std::ostream &os) const
{

}
*/
