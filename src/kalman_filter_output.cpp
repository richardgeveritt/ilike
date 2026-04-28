#include "kalman_filter_output.h"
#include "utils.h"
#include "kalman_filter.h"
#include "filesystem.h"
#include "ilike_hdf5_utils.h"
#include <sstream>

namespace ilike
{
KalmanFilterOutput::KalmanFilterOutput()
:LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  
  this->kf_iteration = 0;
  this->iteration_written_to_file = -1;
  
  this->start_time = std::chrono::high_resolution_clock::now();
}

KalmanFilterOutput::KalmanFilterOutput(KalmanFilter* estimator_in,
                                       size_t lag_in,
                                       const std::string &results_name_in)
:LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  
  if (lag_in<2)
  {
    this->lag = 2;
  }
  else
  {
    this->lag = lag_in;
  }
  this->output_lag = lag_in;
  this->results_name = results_name_in;
  this->kf_iteration = 0;
  this->iteration_written_to_file = -1;
  //this->time = 0.0;
  this->start_time = std::chrono::high_resolution_clock::now();
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
  this->log_normalising_constant_ratios = another.log_normalising_constant_ratios;
  this->schedule_parameters = another.schedule_parameters;
  this->lag = another.lag;
  this->output_lag = another.output_lag;
  
  this->times = another.times;
  this->llhds = another.llhds;
  this->kf_iteration = another.kf_iteration;
  this->iteration_written_to_file = another.iteration_written_to_file;
  this->results_name = another.results_name;
  this->start_time = another.start_time;
  
}

void KalmanFilterOutput::set_time()
{
  size_t num_to_pop_front = std::max<int>(0,this->times.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->times.pop_front();
  }
  
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - this->start_time;
  if (this->times.size()>0)
    this->times.push_back(elapsed_time.count()+this->times.back());
  else
    this->times.push_back(elapsed_time.count());
}

void KalmanFilterOutput::set_llhd(double llhd_in)
{
  size_t num_to_pop_front = std::max<int>(0,this->llhds.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->llhds.pop_front();
  }
  this->llhds.push_back(llhd_in);
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
  size_t num_to_pop_back = std::max<int>(0,this->predicted_means.size()-this->lag+1);
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
  size_t num_to_pop_back = std::max<int>(0,this->posterior_means.size()-this->lag+1);
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

void KalmanFilterOutput::add_schedule_parameters(const Parameters &current_schedule_parameters)
{
  size_t num_to_pop_back = std::max<int>(0,this->schedule_parameters.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->schedule_parameters.pop_back();
  }
  this->schedule_parameters.push_back(current_schedule_parameters);
}

void KalmanFilterOutput::add_log_normalising_constant_ratio(double log_normalising_constant_ratio)
{
  size_t num_to_pop_back = std::max<int>(0,this->log_normalising_constant_ratios.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->log_normalising_constant_ratios.pop_back();
  }
  this->log_normalising_constant_ratios.push_back(log_normalising_constant_ratio);
}

void KalmanFilterOutput::set_current_predicted_to_be_current_posterior()
{
  this->current_predicted_mean = this->current_posterior_mean;
  this->current_predicted_covariance = this->current_posterior_covariance;
}

void KalmanFilterOutput::increment_kf_iteration()
{
  this->kf_iteration = this->kf_iteration + 1;
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
  Rcpp::stop("KalmanFilterOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat KalmanFilterOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                            const Parameters &x)
{
  Rcpp::stop("KalmanFilterOutput::subsample_get_gradient_of_log - not yet implemented.");
}

void KalmanFilterOutput::print(std::ostream &os) const
{
  
}

void KalmanFilterOutput::forget_you_were_already_written_to_file()
{
  this->iteration_written_to_file = -1;
  //int iteration_written_to_file = -1;
}

void KalmanFilterOutput::terminate()
{
  while (this->posterior_means.size()>this->output_lag)
  {
    this->posterior_means.pop_back();
  }
  
  while (this->posterior_covariances.size()>this->output_lag)
  {
    this->posterior_covariances.pop_back();
  }
  
  while (this->llhds.size()>this->output_lag)
  {
    this->llhds.pop_back();
  }
  
  while (this->times.size()>this->output_lag)
  {
    this->times.pop_back();
  }
}

void KalmanFilterOutput::write_to_file(const std::string &dir_name,
                                       const std::string &index)
{
  std::string directory_name = dir_name + "_kf";

  if (!directory_exists(directory_name))
    make_directory(directory_name);

  if (!this->estimator->h5_file)
  {
    this->estimator->h5_file_path = directory_name + "/output.h5";
    this->estimator->h5_file = h5_open_or_create(this->estimator->h5_file_path);

    auto root = this->estimator->h5_file->getGroup("/");
    h5_set_str_attr(root, "variable_names",
                    std::vector<std::string>{this->estimator->state_name});
    std::vector<size_t> vsizes{this->estimator->state_dimension};
    h5_set_sizet_attr(root, "variable_sizes", vsizes);
  }

  HighFive::File &hf = *this->estimator->h5_file;

  for (size_t iteration = this->iteration_written_to_file+1;
       iteration < this->kf_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->kf_iteration - iteration;

    if (this->posterior_means.size() > distance_from_end)
    {
      size_t deque_index = this->posterior_means.size()-1-distance_from_end;

      h5_append_double(hf, "log_likelihood", this->llhds[deque_index]);
      h5_append_double(hf, "time",           this->times[deque_index]);
      h5_append_double(hf, "output_lengths", 1.0);

      std::string iter_grp_path = "iteration/" + std::to_string(iteration+1);
      HighFive::Group iter_grp = h5_ensure_group(hf, iter_grp_path);

      h5_write_scalar_double(iter_grp, "incremental_log_likelihood",
                             this->log_normalising_constant_ratios[deque_index]);
      {
        std::ostringstream oss;
        oss << this->schedule_parameters[deque_index];
        h5_write_string(iter_grp, "schedule_parameters", oss.str());
      }
      h5_write_vec(iter_grp, "posterior_means",
                   this->posterior_means[deque_index]);
      h5_write_mat(iter_grp, "posterior_covariances",
                   this->posterior_covariances[deque_index]);
      h5_write_vec(iter_grp, "predicted_means",
                   this->predicted_means[deque_index]);
      h5_write_mat(iter_grp, "predicted_covariances",
                   this->predicted_covariances[deque_index]);
    }
  }

  this->iteration_written_to_file = this->kf_iteration;
}

void KalmanFilterOutput::close_ofstreams()
{
  this->estimator->h5_file.reset();
}
}
