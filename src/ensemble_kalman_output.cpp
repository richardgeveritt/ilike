#include "ensemble_kalman_output.h"
#include "ensemble_kalman.h"
#include "move_output.h"
#include "filesystem.h"
#include "utils.h"

namespace ilike
{
EnsembleKalmanOutput::EnsembleKalmanOutput()
:LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  this->estimator = NULL;
  this->iteration_written_to_file = -1;
  this->transform = NULL;
  this->enk_iteration = 0;
  this->skip_to_end_of_sequence = false;
  
  this->start_time = std::chrono::high_resolution_clock::now();
}

EnsembleKalmanOutput::EnsembleKalmanOutput(EnsembleKalman* estimator_in,
                                           size_t lag_in,
                                           std::shared_ptr<Transform> transform_in,
                                           const std::string &results_name_in)
:LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in->ensemble_kalman_duplicate();
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
  this->iteration_written_to_file = -1;
  this->results_name = results_name_in;
  this->transform = transform_in;
  this->enk_iteration = 0;
  this->skip_to_end_of_sequence = false;
  this->start_time = std::chrono::high_resolution_clock::now();
}

EnsembleKalmanOutput::~EnsembleKalmanOutput()
{
  if (this->estimator!=NULL)
    delete this->estimator;
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
  
  if (this->estimator!=NULL)
    delete this->estimator;
  
  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* EnsembleKalmanOutput::duplicate() const
{
  return( new EnsembleKalmanOutput(*this));
}

void EnsembleKalmanOutput::make_copy(const EnsembleKalmanOutput &another)
{
  if (another.estimator!=NULL)
    this->estimator = another.estimator->ensemble_kalman_duplicate();
  else
    this->estimator = NULL;
  this->enk_iteration = another.enk_iteration;
  this->iteration_written_to_file = another.iteration_written_to_file;
  this->transform = another.transform;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
  this->all_ensembles = another.all_ensembles;
  this->lag = another.lag;
  this->output_lag = another.output_lag;
  this->start_time = another.start_time;
  this->times = another.times;
  this->llhds = another.llhds;
  this->skip_to_end_of_sequence = another.skip_to_end_of_sequence;
  this->results_name = another.results_name;
}

/*
 Ensemble* EnsembleKalmanOutput::add_ensemble()
 {
 size_t num_to_pop_front = std::max<int>(0,this->all_ensembles.size()-this->lag+1);
 for (size_t i=0; i<num_to_pop_front; ++i)
 {
 this->all_ensembles.pop_front();
 }
 this->all_ensembles.push_back(Ensemble());
 this->all_ensembles.back().reserve(this->estimator->number_of_ensemble_members);
 return &this->all_ensembles.back();
 }
 */

Ensemble* EnsembleKalmanOutput::add_ensemble(EnsembleFactors* ensemble_factors)
{
  size_t num_to_pop_front = std::max<int>(0,this->all_ensembles.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->all_ensembles.pop_front();
  }
  this->all_ensembles.push_back(Ensemble(ensemble_factors));
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

double EnsembleKalmanOutput::calculate_latest_log_normalising_constant_ratio()
{
  return this->all_ensembles.back().calculate_log_normalising_constant();
}

double EnsembleKalmanOutput::calculate_inversion_latest_log_normalising_constant_ratio(double inverse_incremental_temperature)
{
  return this->all_ensembles.back().calculate_inversion_log_normalising_constant(inverse_incremental_temperature);
}

double EnsembleKalmanOutput::calculate_unbiased_inversion_latest_log_normalising_constant_ratio(double inverse_incremental_temperature)
{
  return this->all_ensembles.back().calculate_unbiased_inversion_log_normalising_constant(inverse_incremental_temperature);
}

double EnsembleKalmanOutput::calculate_path1_inversion_latest_log_normalising_constant_ratio(const std::vector<double> &previous_log_measurement_likelihood_means,
                                                                                             double inverse_incremental_temperature,
                                                                                             double temperature,
                                                                                             double multiplier)
{
  this->all_ensembles.back().calculate_path1_inversion_log_normalising_constant((this->all_ensembles.end()-1)->log_measurement_likelihood_means,
                                                                                temperature,
                                                                                multiplier);
  
  std::deque<Ensemble>::iterator current_ensemble = this->all_ensembles.end()-1;
  //std::deque<Ensemble>::iterator previous_ensemble = this->all_ensembles.end()-2;
  
  double log_ratio = 0.0;
  for (size_t i=0; i<current_ensemble->log_measurement_likelihood_means.size(); ++i)
  {
    log_ratio = log_ratio + (0.5/inverse_incremental_temperature)*(current_ensemble->log_measurement_likelihood_means[i] + previous_log_measurement_likelihood_means[i]);
  }
  
  return log_ratio;
}

double EnsembleKalmanOutput::calculate_path2_inversion_latest_log_normalising_constant_ratio(const std::vector<double> &previous_log_measurement_likelihood_means,
                                                                                             const std::vector<double> &previous_log_measurement_likelihood_variances,
                                                                                             double inverse_incremental_temperature)
{
  this->all_ensembles.back().calculate_path2_inversion_log_normalising_constant((this->all_ensembles.end()-1)->log_measurement_likelihood_means,
                                                                                (this->all_ensembles.end()-1)->log_measurement_likelihood_variances);
  
  std::deque<Ensemble>::iterator current_ensemble = this->all_ensembles.end()-1;
  //std::deque<Ensemble>::iterator previous_ensemble = this->all_ensembles.end()-2;
  
  double log_ratio = 0.0;
  for (size_t i=0; i<current_ensemble->log_measurement_likelihood_means.size(); ++i)
  {
    log_ratio = log_ratio + (1.0/(2.0*inverse_incremental_temperature))*(current_ensemble->log_measurement_likelihood_means[i] + previous_log_measurement_likelihood_means[i]);
    
    log_ratio = log_ratio - (1.0/(12.0*pow(inverse_incremental_temperature,2.0)))*(current_ensemble->log_measurement_likelihood_variances[i] - previous_log_measurement_likelihood_variances[i]);
  }
  
  return log_ratio;
}

void EnsembleKalmanOutput::calculate_path1_inversion_initial_quantities(std::vector<double> &initial_log_measurement_likelihood_means,
                                                                        double temperature,
                                                                        double multiplier)
{
  this->all_ensembles.back().calculate_path1_inversion_log_normalising_constant(initial_log_measurement_likelihood_means,
                                                                                temperature,
                                                                                multiplier);
}

void EnsembleKalmanOutput::calculate_path2_inversion_initial_quantities(std::vector<double> &initial_log_measurement_likelihood_means,
                                                                        std::vector<double> &initial_log_measurement_likelihood_variances)
{
  this->all_ensembles.back().calculate_path2_inversion_log_normalising_constant(initial_log_measurement_likelihood_means,
                                                                                initial_log_measurement_likelihood_variances);
}

void EnsembleKalmanOutput::calculate_kalman_gains(double inverse_incremental_temperature)
{
  this->all_ensembles.back().calculate_kalman_gains(inverse_incremental_temperature);
}

void EnsembleKalmanOutput::simulate()
{
  this->estimator->ensemble_kalman_simulate(this);
}

void EnsembleKalmanOutput::simulate(const Parameters &parameters)
{
  this->estimator->ensemble_kalman_simulate(this, parameters);
}

/*
 void EnsembleKalmanOutput::evaluate_smcfixed_part()
 {
 if (this->estimator->smcfixed_flag)
 {
 this->estimator->ensemble_kalman_evaluate(this);
 //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
 //this->log_normalising_constant_ratios.end(),
 //1.0);
 }
 else
 {
 this->estimator->ensemble_kalman_evaluate_smcfixed_part(this);
 }
 }
 
 void EnsembleKalmanOutput::evaluate_smcadaptive_part_given_smcfixed()
 {
 if (!this->estimator->smcfixed_flag)
 {
 this->estimator->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(this);
 }
 
 }
 */

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

void EnsembleKalmanOutput::subsample_simulate()
{
  this->estimator->ensemble_kalman_simulate(this);
}

void EnsembleKalmanOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator->ensemble_kalman_simulate(this, parameters);
}

/*
 void EnsembleKalmanOutput::subsample_evaluate_smcfixed_part()
 {
 if (this->estimator->smcfixed_flag)
 {
 this->estimator->ensemble_kalman_subsample_evaluate(this);
 //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
 //this->log_normalising_constant_ratios.end(),
 //1.0);
 }
 else
 {
 this->estimator->ensemble_kalman_subsample_evaluate_smcfixed_part(this);
 }
 }
 
 void EnsembleKalmanOutput::subsample_evaluate_smcadaptive_part_given_smcfixed()
 {
 if (!this->estimator->smcfixed_flag)
 {
 this->estimator->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(this);
 }
 
 }
 */

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
  Rcpp::stop("EnsembleKalmanOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat EnsembleKalmanOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                              const Parameters &x)
{
  Rcpp::stop("EnsembleKalmanOutput::subsample_get_gradient_of_log - not yet implemented.");
}

void EnsembleKalmanOutput::set_time()
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

void EnsembleKalmanOutput::set_llhd(double llhd_in)
{
  size_t num_to_pop_front = std::max<int>(0,this->llhds.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->llhds.pop_front();
  }
  this->llhds.push_back(llhd_in);
}

size_t EnsembleKalmanOutput::number_of_ensemble_kalman_iterations() const
{
  return this->all_ensembles.size();
}

void EnsembleKalmanOutput::increment_enk_iteration()
{
  this->enk_iteration = this->enk_iteration + 1;
}

void EnsembleKalmanOutput::forget_you_were_already_written_to_file()
{
  this->iteration_written_to_file = -1;
}

void EnsembleKalmanOutput::terminate()
{
  while (this->all_ensembles.size()>this->output_lag)
  {
    this->all_ensembles.pop_back();
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

void EnsembleKalmanOutput::skip_to_end_of_sequence_if_points_are_gaussian(double significance_level)
{
  if (hz(this->back().get_packed_members())>significance_level)
  {
    this->skip_to_end_of_sequence = true;
  }
}

void EnsembleKalmanOutput::write_to_file(const std::string &dir_name,
                                         const std::string &index)
{
  std::string directory_name = dir_name + "_enk";
  
  //if (index!="")
  //  directory_name = directory_name + "_" + index;
  
  if (!directory_exists(directory_name))
  {
    make_directory(directory_name);
  }
  
  // for each iteration left to write
  for (size_t iteration = this->iteration_written_to_file+1;
       iteration<this->enk_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->enk_iteration-iteration;
    
    //size_t llhd_index = this->llhds.size()-1-distance_from_end;
    
    if (this->all_ensembles.size() > distance_from_end)
    {
      size_t deque_index = this->all_ensembles.size()-1-distance_from_end;
      
      if (!this->estimator->log_likelihood_file_stream.is_open())
      {
        this->estimator->log_likelihood_file_stream.open(directory_name + "/log_likelihood.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->log_likelihood_file_stream.is_open())
      {
        this->estimator->log_likelihood_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->llhds[deque_index] << std::endl;
        //log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/log_likelihood.txt" + "cannot be opened.");
      }
      
      if (!this->estimator->time_file_stream.is_open())
      {
        this->estimator->time_file_stream.open(directory_name + "/time.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->time_file_stream.is_open())
      {
        //double time_sum = 0.0;
        //for (size_t k=0; k<=llhd_index; ++k)
        //{
        //  time_sum = time_sum + this->times[k];
        //}
        this->estimator->time_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->times[deque_index] << std::endl;
        //log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/time.txt" + " cannot be opened.");
      }
      
      
      if (!this->estimator->vector_variables_file_stream.is_open())
      {
        this->estimator->vector_variables_file_stream.open(directory_name + "/vector_variables.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->vector_variables_file_stream.is_open())
      {
        for (size_t i=0; i<this->estimator->vector_variables.size(); ++i)
        {
          this->estimator->vector_variables_file_stream << this->estimator->vector_variables[i] << ";";
        }
        this->estimator->vector_variables_file_stream << std::endl;
        //vector_variables_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/vector_variables.txt" + "cannot be opened.");
      }
      
      if (!this->estimator->vector_variable_sizes_file_stream.is_open())
      {
        this->estimator->vector_variable_sizes_file_stream.open(directory_name + "/vector_variable_sizes.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->vector_variable_sizes_file_stream.is_open())
      {
        for (size_t i=0; i<this->estimator->vector_variable_sizes.size(); ++i)
        {
          this->estimator->vector_variable_sizes_file_stream << this->estimator->vector_variable_sizes[i] << ";";
        }
        this->estimator->vector_variable_sizes_file_stream << std::endl;
        //vector_variable_sizes_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/vector_variable_sizes.txt" + "cannot be opened.");
      }
      
      if (!this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream.open(directory_name + "/output_lengths.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream << this->all_ensembles[deque_index].get_output_lengths();
        //incremental_log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + directory_name + "/output_lengths.txt" + " cannot be opened.");
      }
      
      std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(iteration+1);
      
      if (!directory_exists(smc_iteration_directory))
      {
        make_directory(smc_iteration_directory);
      }
      
      if (!this->estimator->incremental_log_likelihood_file_stream.is_open())
      {
        this->estimator->incremental_log_likelihood_file_stream.open(smc_iteration_directory + "/incremental_log_likelihood.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->incremental_log_likelihood_file_stream.is_open())
      {
        this->estimator->incremental_log_likelihood_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->all_ensembles[deque_index].log_normalising_constant_ratio << std::endl;
        //incremental_log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/incremental_log_likelihood.txt" + "cannot be opened.");
      }
      
      if (!this->estimator->ess_file_stream.is_open())
      {
        this->estimator->ess_file_stream.open(smc_iteration_directory + "/ess.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->ess_file_stream.is_open())
      {
        this->estimator->ess_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->all_ensembles[deque_index].ess << std::endl;
        //ess_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/ess.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->schedule_parameters_file_stream.is_open())
      {
        this->estimator->schedule_parameters_file_stream.open(smc_iteration_directory + "/schedule_parameters.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->schedule_parameters_file_stream.is_open())
      {
        if (this->estimator->reciprocal_schedule_scale==0.0)
        {
          this->estimator->schedule_parameters_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->all_ensembles[deque_index].schedule_parameters << std::endl;
        }
        else
        {
          this->estimator->schedule_parameters_file_stream << std::setprecision(std::numeric_limits<double>::max_digits10) << this->estimator->reciprocal_schedule_scale * pow(this->all_ensembles[deque_index].schedule_parameters,-0.5) << std::endl;
        }
        //schedule_parameters_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/schedule_parameters.txt" + "cannot be opened.");
      }
      
      if(!this->estimator->vector_points_file_stream.is_open())
      {
        this->estimator->vector_points_file_stream.open(smc_iteration_directory + "/vector_points.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->vector_points_file_stream.is_open())
      {
        for (auto i = this->all_ensembles[deque_index].members.begin();
             i!=this->all_ensembles[deque_index].members.end();
             ++i)
        {
          (*i)->write_vector_points(this->estimator->vector_variables,
                                    this->estimator->vector_points_file_stream,
                                    this->transform);
        }
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/vector_points.txt" + "cannot be opened.");
      }
      
      if(!this->estimator->any_points_file_stream.is_open())
      {
        this->estimator->any_points_file_stream.open(smc_iteration_directory + "/any_points.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->any_points_file_stream.is_open())
      {
        for (auto i = this->all_ensembles[deque_index].members.begin();
             i!=this->all_ensembles[deque_index].members.end();
             ++i)
        {
          (*i)->write_any_points(this->estimator->any_variables,
                                 this->estimator->any_points_file_stream);
        }
        //any_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/any_points.txt" + "cannot be opened.");
      }
      
      // for any other info in Particles, write to a different file
      
      //for (auto i = this->all_ensembles[deque_index].members.begin();
      //     i!=this->all_ensembles[deque_index].members.end();
      //     ++i)
      //{
      //  (*i)->write_factors(smc_iteration_directory);
      //}
      for (size_t i = 0;
           i<this->all_ensembles[deque_index].members.size();
           ++i)
      {
        this->all_ensembles[deque_index].members[i]->write_factors(smc_iteration_directory,
                                                                   std::to_string(i));
      }
      
      //for (auto i = this->all_ensembles[deque_index].members.begin();
      //     i!=this->all_ensembles[deque_index].members.end();
      //     ++i)
      //{
      //  (*i)->write_ensemble_factors(smc_iteration_directory);
      //}
      for (size_t i = 0;
           i<this->all_ensembles[deque_index].members.size();
           ++i)
      {
        this->all_ensembles[deque_index].members[i]->write_ensemble_factors(smc_iteration_directory,
                                                                            std::to_string(i));
      }
      
      this->close_ofstreams(deque_index);
      
    }
    
  }
  
  this->iteration_written_to_file = this->enk_iteration;
}

void EnsembleKalmanOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->time_file_stream.close();
  this->estimator->vector_variables_file_stream.close();
  this->estimator->vector_variable_sizes_file_stream.close();
  this->estimator->incremental_log_likelihood_file_stream.close();
  this->estimator->output_lengths_file_stream.close();
  this->estimator->ess_file_stream.close();
  this->estimator->schedule_parameters_file_stream.close();
  this->estimator->vector_points_file_stream.close();
  this->estimator->any_points_file_stream.close(); // should be one for each member of Parameters
  this->estimator->time_file_stream.close();
  
  for (auto i = this->all_ensembles.begin();
       i!=this->all_ensembles.end();
       ++i)
  {
    i->close_ofstreams();
  }
}

void EnsembleKalmanOutput::close_ofstreams(size_t deque_index)
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->time_file_stream.close();
  this->estimator->vector_variables_file_stream.close();
  this->estimator->vector_variable_sizes_file_stream.close();
  this->estimator->incremental_log_likelihood_file_stream.close();
  this->estimator->output_lengths_file_stream.close();
  this->estimator->ess_file_stream.close();
  this->estimator->schedule_parameters_file_stream.close();
  this->estimator->vector_points_file_stream.close();
  this->estimator->any_points_file_stream.close(); // should be one for each member of Parameters
  this->estimator->time_file_stream.close();
  
  this->all_ensembles[deque_index].close_ofstreams();
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
 size_t num_to_pop_front = std::max<int>(0,this->predicted_means.size()-this->lag+1);
 for (size_t i=0; i<num_to_pop_front; ++i)
 {
 this->predicted_means.pop_front();
 }
 this->predicted_means.push_back(this->current_predicted_mean);
 for (size_t i=0; i<num_to_pop_front; ++i)
 {
 this->predicted_covariances.pop_front();
 }
 this->predicted_covariances.push_back(this->current_predicted_covariance);
 }
 
 void EnsembleKalmanOutput::add_posterior_statistics()
 {
 // unsure
 // copied from KF
 size_t num_to_pop_front = std::max<int>(0,this->posterior_means.size()-this->lag+1);
 for (size_t i=0; i<num_to_pop_front; ++i)
 {
 this->posterior_means.pop_front();
 }
 this->posterior_means.push_back(this->current_posterior_mean);
 for (size_t i=0; i<num_to_pop_front; ++i)
 {
 this->posterior_covariances.pop_front();
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
}
