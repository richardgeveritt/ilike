#include "kalman_filter_output.h"
#include "utils.h"
#include "kalman_filter.h"

KalmanFilterOutput::KalmanFilterOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

KalmanFilterOutput::KalmanFilterOutput(KalmanFilter* estimator_in,
                                       size_t lag_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  
  this->lag = lag_in;
}

KalmanFilterOutput::~KalmanFilterOutput()
{

}

//Copy constructor for the KalmanFilterOutput class.
KalmanFilterOutput::KalmanFilterOutput(const KalmanFilterOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void KalmanFilterOutput::operator=(const KalmanFilterOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* KalmanFilterOutput::duplicate() const
{
  return( new KalmanFilterOutput(*this));
}

void KalmanFilterOutput::make_copy(const KalmanFilterOutput &another)
{
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->posterior_means = another.posterior_means;
  this->posterior_covariances = another.posterior_covariances;
  this->predicted_means = another.predicted_means;
  this->predicted_covariances = another.predicted_covariances;
  this->current_predicted_mean = another.current_predicted_mean;
  this->current_predicted_covariance = another.current_predicted_covariance;
  this->current_posterior_mean = another.current_posterior_mean;
  this->current_posterior_covariance = another.current_posterior_covariance;
  this->lag = another.lag;
}

void KalmanFilterOutput::simulate()
{
  // Deterministic, so nothing happens here.
}


void KalmanFilterOutput::simulate(const Parameters &parameters)
{
  // Deterministic, so nothing happens here.
}

void KalmanFilterOutput::evaluate_smcfixed_part(const Parameters &conditioned_on_parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate(this, conditioned_on_parameters);
    this->log_likelihood_smcfixed_part = this->log_likelihood;
  }
}

void KalmanFilterOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate(this, conditioned_on_parameters);
  }
  else
  {
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  }
  
}

void KalmanFilterOutput::subsample_simulate(const Parameters &parameters)
{
  // Deterministic, so nothing happens here.
}

void KalmanFilterOutput::subsample_evaluate_smcfixed_part(const Parameters &conditioned_on_parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate(this,conditioned_on_parameters);
    this->subsample_log_likelihood_smcfixed_part = this->subsample_log_likelihood;
  }
}

void KalmanFilterOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate(this, conditioned_on_parameters);
  }
  else
  {
    this->subsample_log_likelihood = this->subsample_log_likelihood_smcfixed_part;
  }
  
}

void KalmanFilterOutput::set_current_predicted_statistics(const arma::colvec &latest_mean,
                                                          const arma::mat &latest_covariance)
{
  this->current_predicted_mean = latest_mean;
  this->current_predicted_covariance = latest_covariance;
}

void KalmanFilterOutput::set_current_posterior_statistics(const arma::colvec &latest_mean,
                                                          const arma::mat &latest_covariance)
{
  this->current_posterior_mean = latest_mean;
  this->current_posterior_covariance = latest_covariance;
}

void KalmanFilterOutput::add_predicted_statistics()
{
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

void KalmanFilterOutput::add_posterior_statistics()
{
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

void KalmanFilterOutput::set_current_predicted_to_be_current_posterior()
{
  this->current_predicted_mean = this->current_posterior_mean;
  this->current_predicted_covariance = this->current_posterior_covariance;
}

LikelihoodEstimator* KalmanFilterOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::colvec KalmanFilterOutput::predicted_mean_back() const
{
  return this->predicted_means.back();
}

arma::colvec KalmanFilterOutput::posterior_mean_back() const
{
  return this->posterior_means.back();
}

arma::mat KalmanFilterOutput::predicted_covariance_back() const
{
  return this->predicted_covariances.back();
}

arma::mat KalmanFilterOutput::posterior_covariance_back() const
{
  return this->posterior_covariances.back();
}

size_t KalmanFilterOutput::predicted_size() const
{
  return this->predicted_means.size();
}

arma::mat KalmanFilterOutput::get_gradient_of_log(const std::string &variable,
                                                  const Parameters &x)
{
  throw std::runtime_error("KalmanFilterOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat KalmanFilterOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                  const Parameters &x)
{
  throw std::runtime_error("KalmanFilterOutput::subsample_get_gradient_of_log - not yet implemented.");
}

void KalmanFilterOutput::print(std::ostream &os) const
{

}
