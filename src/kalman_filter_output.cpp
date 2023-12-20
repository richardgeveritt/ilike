#include "kalman_filter_output.h"
#include "utils.h"
#include "kalman_filter.h"
#include "filesystem.h"

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
  
  this->lag = lag_in;
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
  this->schedule_parameters = another.schedule_parameters;
  this->lag = another.lag;
  
  this->times = another.times;
  this->llhds = another.llhds;
  this->kf_iteration = another.kf_iteration;
  this->iteration_written_to_file = another.iteration_written_to_file;
  this->results_name = another.results_name;
  this->start_time = another.start_time;
  
}

void KalmanFilterOutput::set_time()
{
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - this->start_time;
  this->times.push_back(elapsed_time.count());
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

void KalmanFilterOutput::write_to_file(const std::string &dir_name,
                                       const std::string &index)
{
  std::string directory_name = dir_name + "_kf";
  
  //if (index!="")
  //  directory_name = directory_name + "_" + index;
  
  if (!directory_exists(directory_name))
  {
    make_directory(directory_name);
  }
  
  // for each iteration left to write
  for (size_t iteration = this->iteration_written_to_file+1;
       iteration<this->kf_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->kf_iteration-iteration;
    
    size_t llhd_index = this->llhds.size()-1-distance_from_end;
    
    //if (int(this->llhds.size())-1-int(distance_from_end)>=0)
    //{
    if (!this->estimator->log_likelihood_file_stream.is_open())
    {
      this->estimator->log_likelihood_file_stream.open(directory_name + "/log_likelihood.txt",std::ios::out | std::ios::app);
    }
    if (this->estimator->log_likelihood_file_stream.is_open())
    {
      this->estimator->log_likelihood_file_stream << this->llhds[llhd_index] << std::endl;
      //log_likelihood_file_stream.close();
    }
    else
    {
      Rcpp::stop("File " + directory_name + "/log_likelihood.txt" + " cannot be opened.");
    }
    //}
    
    if (!this->estimator->time_file_stream.is_open())
    {
      this->estimator->time_file_stream.open(directory_name + "/time.txt",std::ios::out | std::ios::app);
    }
    if (this->estimator->time_file_stream.is_open())
    {
      double time_sum = 0.0;
      for (size_t k=0; k<this->times.size(); ++k)
      {
        time_sum = time_sum + this->times[k];
      }
      this->estimator->time_file_stream << time_sum << std::endl;
      //log_likelihood_file_stream.close();
    }
    else
    {
      Rcpp::stop("File " + directory_name + "/time.txt" + " cannot be opened.");
    }
    
    if (this->posterior_means.size() > distance_from_end)
    {
      size_t deque_index = this->posterior_means.size()-1-distance_from_end;
      
      if (!this->estimator->vector_variables_file_stream.is_open())
      {
        this->estimator->vector_variables_file_stream.open(directory_name + "/vector_variables.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->vector_variables_file_stream.is_open())
      {
        this->estimator->vector_variables_file_stream << this->estimator->state_name << ";";
        this->estimator->vector_variables_file_stream << std::endl;
        //vector_variables_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/vector_variables.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->vector_variable_sizes_file_stream.is_open())
      {
        this->estimator->vector_variable_sizes_file_stream.open(directory_name + "/vector_variable_sizes.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->vector_variable_sizes_file_stream.is_open())
      {
        this->estimator->vector_variable_sizes_file_stream << this->estimator->state_dimension << ";";
        this->estimator->vector_variable_sizes_file_stream << std::endl;
        //vector_variable_sizes_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/vector_variable_sizes.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream.open(directory_name + "/output_lengths.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream << 1 << std::endl;
        //incremental_log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/output_lengths.txt" + " cannot be opened.");
      }
      
      /*
       std::ofstream any_variables_file_stream;
       any_variables_file_stream.open(directory_name + "/any_variables.txt",std::ios::out | std::ios::app);
       if (any_variables_file_stream.is_open())
       {
       for (size_t i=0; i<this->estimator->any_variables.size(); ++i)
       {
       any_variables_file_stream << this->estimator->any_variables[i] << ";";
       }
       any_variables_file_stream << std::endl;
       any_variables_file_stream.close();
       }
       else
       {
       Rcpp::stop("File " + directory_name + "/any_variables.txt" + " cannot be opened.");
       }
       */
      
      std::string kf_iteration_directory = directory_name + "/iteration" + std::to_string(iteration+1);
      
      if (!directory_exists(kf_iteration_directory))
      {
        make_directory(kf_iteration_directory);
      }
      
      if (!this->estimator->incremental_log_likelihood_file_stream.is_open())
      {
        this->estimator->incremental_log_likelihood_file_stream.open(kf_iteration_directory + "/incremental_log_likelihood.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->incremental_log_likelihood_file_stream.is_open())
      {
        this->estimator->incremental_log_likelihood_file_stream << this->log_normalising_constant_ratios[deque_index] << std::endl;
        //incremental_log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/incremental_log_likelihood.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->schedule_parameters_file_stream.is_open())
      {
        this->estimator->schedule_parameters_file_stream.open(kf_iteration_directory + "/schedule_parameters.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->schedule_parameters_file_stream.is_open())
      {
        this->estimator->schedule_parameters_file_stream << this->schedule_parameters[deque_index] << std::endl;
        //schedule_parameters_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/schedule_parameters.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->posterior_means_file_stream.is_open())
      {
        this->estimator->posterior_means_file_stream.open(kf_iteration_directory + "/posterior_means.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->posterior_means_file_stream.is_open())
      {
        this->estimator->posterior_means_file_stream << this->posterior_means[deque_index].as_row();
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/posterior_means.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->posterior_covariances_file_stream.is_open())
      {
        this->estimator->posterior_covariances_file_stream.open(kf_iteration_directory + "/posterior_covariances.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->posterior_covariances_file_stream.is_open())
      {
        this->estimator->posterior_covariances_file_stream << this->posterior_covariances[deque_index].as_row();
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/posterior_covariances.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->predicted_means_file_stream.is_open())
      {
        this->estimator->predicted_means_file_stream.open(kf_iteration_directory + "/predicted_means.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->predicted_means_file_stream.is_open())
      {
        this->estimator->predicted_means_file_stream << this->predicted_means[deque_index].as_row();
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/predicted_means.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->predicted_covariances_file_stream.is_open())
      {
        this->estimator->predicted_covariances_file_stream.open(kf_iteration_directory + "/predicted_covariances.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->predicted_covariances_file_stream.is_open())
      {
        this->estimator->predicted_covariances_file_stream << this->predicted_covariances[deque_index].as_row();
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + kf_iteration_directory + "/predicted_covariances.txt" + " cannot be opened.");
      }
      
      //this->close_ofstreams(deque_index);
      this->close_ofstreams();
    }
  }
  
  this->iteration_written_to_file = this->kf_iteration;
}

void KalmanFilterOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->time_file_stream.close();
  this->estimator->vector_variables_file_stream.close();
  this->estimator->vector_variable_sizes_file_stream.close();
  this->estimator->output_lengths_file_stream.close();
  this->estimator->incremental_log_likelihood_file_stream.close();
  this->estimator->schedule_parameters_file_stream.close();
  this->estimator->posterior_means_file_stream.close();
  this->estimator->posterior_covariances_file_stream.close();
  this->estimator->predicted_means_file_stream.close();
  this->estimator->predicted_covariances_file_stream.close();
}
