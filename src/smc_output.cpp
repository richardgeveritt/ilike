#include <algorithm>
#include <numeric>
#include "smc_output.h"
#include "smc.h"
#include "importance_sampler.h"
#include "smc_mcmc_move.h"
#include "smc_marginal.h"
#include "smc_generic.h"
#include "filesystem.h"
#include "move_output.h"

SMCOutput::SMCOutput()
  :LikelihoodEstimatorOutput()
{
  this->estimator = NULL;
  this->smc_iteration = 0;
  this->iteration_written_to_file = -1;
  //this->time = 0.0;
  
  this->start_time = std::chrono::high_resolution_clock::now();
}

SMCOutput::~SMCOutput()
{

}

SMCOutput::SMCOutput(SMC* estimator_in,
                     size_t lag_in,
                     size_t lag_proposed_in,
                     const std::string &results_name_in)
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_pre_last_step = 0.0;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->estimator = estimator_in;
  this->results_name = results_name_in;
  this->smc_iteration = 0;
  this->iteration_written_to_file = -1;
  //this->time = 0.0;
  this->start_time = std::chrono::high_resolution_clock::now();
}

//Copy constructor for the SMCOutput class.
SMCOutput::SMCOutput(const SMCOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void SMCOutput::operator=(const SMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->all_particles.clear();
  this->all_proposed.clear();
  //this->unnormalised_log_weights.clear();
  //this->normalised_log_weights.clear();
  //this->log_normalising_constant_ratios.clear();
  //this->incremental_log_weights.clear();

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* SMCOutput::duplicate() const
{
  return( new SMCOutput(*this));
}

SMCOutput* SMCOutput::smc_duplicate(void)const
{
  return( new SMCOutput(*this));
}

void SMCOutput::make_copy(const SMCOutput &another)
{
  this->all_particles = another.all_particles;
  this->all_proposed = another.all_proposed;
  //this->log_normalising_constant_ratios = another.log_normalising_constant_ratios;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  this->estimator = another.estimator;
  this->log_likelihood_pre_last_step = another.log_likelihood_pre_last_step;
  this->results_name = another.results_name;
  this->smc_iteration = another.smc_iteration;
  this->iteration_written_to_file = another.iteration_written_to_file;
  this->start_time = another.start_time;
  this->times = another.times;
  this->llhds = another.llhds;
}

void SMCOutput::simulate()
{
  this->estimator->simulate_smc(this);
}

void SMCOutput::evaluate_smcfixed_part()
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smc(this);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->evaluate_smcfixed_part_smc(this);
  }
}

void SMCOutput::evaluate_smcadaptive_part_given_smcfixed()
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smcadaptive_part_given_smcfixed_smc(this);
  }
}

void SMCOutput::simulate(const Parameters &parameters)
{
  this->estimator->simulate_smc(this, parameters);
}

void SMCOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smc(this, parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
                        //this->log_normalising_constant_ratios.end(),
                        //1.0);
  }
  else
  {
    this->estimator->evaluate_smcfixed_part_smc(this, parameters);
  }
}

void SMCOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smcadaptive_part_given_smcfixed_smc(this,parameters);
  }
}

void SMCOutput::subsample_simulate()
{
  this->estimator->subsample_simulate_smc(this);
}

void SMCOutput::subsample_evaluate_smcfixed_part()
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smc(this);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->subsample_evaluate_smcfixed_part_smc(this);
  }
}

void SMCOutput::subsample_evaluate_smcadaptive_part_given_smcfixed()
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(this);
  }
}

void SMCOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator->subsample_simulate_smc(this, parameters);
}

void SMCOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smc(this, parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->subsample_evaluate_smcfixed_part_smc(this, parameters);
  }
}

void SMCOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(this,parameters);
  }
}

Particles* SMCOutput::add_particles()
{
  size_t num_to_pop_front = std::max<int>(0,this->all_particles.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->all_particles.pop_front();
  }
  this->all_particles.push_back(Particles());
  this->all_particles.back().reserve(this->estimator->number_of_particles);
  return &this->all_particles.back();
}

Particles* SMCOutput::add_particles(Particles* most_recent_particles)
{
  size_t num_to_pop_front = std::max<int>(0,this->all_particles.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->all_particles.pop_front();
  }
  this->all_particles.push_back(Particles());
  this->all_particles.back().reserve(this->estimator->number_of_particles);
  
  if ((this->all_particles.end()-2)->resampled_flag==true)
  {
    this->all_particles.back().previous_normalised_log_weights.fill(-log(double(this->estimator->number_of_particles)));
  }
  else
  {
    this->all_particles.back().previous_normalised_log_weights = most_recent_particles->normalised_log_weights;
  }
  
  return &this->all_particles.back();
}

void SMCOutput::add_proposed_particles(const Particles &latest_proposed_particles)
{
  size_t num_to_pop_front = std::max<int>(0,this->all_proposed.size()-lag_proposed-1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->all_proposed.pop_front();
  }
  this->all_proposed.push_back(latest_proposed_particles);
}

Particles SMCOutput::back() const
{
  return this->all_particles.back();
}

Particles& SMCOutput::back()
{
  return this->all_particles.back();
}

std::deque<Particles>::iterator SMCOutput::end()
{
  return this->all_particles.end();
}

std::deque<Particles>::const_iterator SMCOutput::end() const
{
  return this->all_particles.end();
}

double SMCOutput::calculate_latest_log_normalising_constant_ratio()
{
  return this->all_particles.back().calculate_log_normalising_constant();
}

void SMCOutput::update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
{
  this->all_particles.back().update_weights(latest_unnormalised_log_incremental_weights);
}

//void SMCOutput::initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
//{
  //arma::colvec latest_unnormalised_log_weights;
  //if (this->unnormalised_log_incremental_weights.size()>0)
  //  latest_unnormalised_log_weights = this->unnormalised_log_weights.back() + latest_unnormalised_log_incremental_weights;
  //else
  //  latest_unnormalised_log_weights = latest_unnormalised_log_weight_updates;
  //this->unnormalised_log_incremental_weights.push_back(latest_unnormalised_log_incremental_weights);
  //size_t num_to_pop_front = std::max<int>(0,unnormalised_log_incremental_weights.size()-lag);
  //for (size_t i=0; i<num_to_pop_front; ++i)
  //{
  //  this->unnormalised_log_incremental_weights.pop_front();
  //}
//}

//void SMCOutput::set_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
//{
//  this->unnormalised_log_incremental_weights.push_back(latest_unnormalised_log_incremental_weights);
//  size_t num_to_pop_front = std::max<int>(0,this->unnormalised_log_incremental_weights.size()-lag);
//  for (size_t i=0; i<num_to_pop_front; ++i)
//  {
//    this->unnormalised_log_incremental_weights.pop_front();
//  }
//}

//void SMCOutput::initialise_next_step()
//{
//  arma::colvec init(this->estimator->number_of_particles);
//  init.fill(0.0);
//  this->unnormalised_log_weights.push_back(init);
//  this->incremental_log_weights.push_back(init);
//  //this->unnormalised_log_incremental_weights.push_back(init);
//  this->log_normalising_constant_ratios.push_back(0.0);
//}

void SMCOutput::normalise_and_resample_weights()
{
  this->log_likelihood_pre_last_step = this->log_likelihood;
  this->all_particles.back().normalise_weights();
  this->resample();
  this->set_time();
  
  if (this->results_name!="")
    this->write(results_name);
  
  this->start_time = std::chrono::high_resolution_clock::now();
}

void SMCOutput::resample()
{
  this->estimator->resample(this);
}

LikelihoodEstimator* SMCOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

size_t SMCOutput::number_of_smc_iterations() const
{
  return this->all_particles.size();
}

arma::mat SMCOutput::get_gradient_of_log(const std::string &variable,
                                         const Parameters &x)
{
  Rcpp::stop("SMCOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat SMCOutput::subsample_get_gradient_of_log(const std::string &variable,
                                         const Parameters &x)
{
  Rcpp::stop("SMCOutput::get_gradient_of_log - not yet implemented.");
}

void SMCOutput::set_time()
{
  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - this->start_time;
  this->times.push_back(elapsed_time.count());
}

void SMCOutput::forget_you_were_already_written_to_file()
{
  this->iteration_written_to_file = -1;
}

void SMCOutput::write_to_file(const std::string &dir_name,
                              const std::string &index)
{
  std::string directory_name = dir_name + "_smc";
  
  //if (index!="")
  //  directory_name = directory_name + "_" + index;
  
  if (!directory_exists(directory_name))
  {
    make_directory(directory_name);
  }
  
  // for each iteration left to write
  for (size_t iteration = this->iteration_written_to_file+1;
       iteration<this->smc_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->smc_iteration-iteration;
   
    size_t llhd_index = this->llhds.size()-1-distance_from_end;
    
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
    
    if (this->all_particles.size() > distance_from_end)
    {
      size_t deque_index = this->all_particles.size()-1-distance_from_end;
      
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
        Rcpp::stop("File " + directory_name + "/vector_variables.txt" + " cannot be opened.");
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
        Rcpp::stop("File " + directory_name + "/vector_variable_sizes.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream.open(directory_name + "/output_lengths.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->output_lengths_file_stream.is_open())
      {
        this->estimator->output_lengths_file_stream << this->all_particles[deque_index].get_output_lengths();
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
        this->estimator->incremental_log_likelihood_file_stream << this->all_particles[deque_index].log_normalising_constant_ratio << std::endl;
        //incremental_log_likelihood_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/incremental_log_likelihood.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->resampled_file_stream.is_open())
      {
        this->estimator->resampled_file_stream.open(smc_iteration_directory + "/resampled.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->resampled_file_stream.is_open())
      {
        this->estimator->resampled_file_stream << this->all_particles[deque_index].resampled_flag << std::endl;
        //resampled_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/resampled.txt" + " cannot be opened.");
      }
      
      if (!this->estimator->ess_file_stream.is_open())
      {
        this->estimator->ess_file_stream.open(smc_iteration_directory + "/ess.txt",std::ios::out | std::ios::app);
      }
      if (this->estimator->ess_file_stream.is_open())
      {
        this->estimator->ess_file_stream << this->all_particles[deque_index].ess << std::endl;
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
        this->estimator->schedule_parameters_file_stream << this->all_particles[deque_index].schedule_parameters << std::endl;
        //schedule_parameters_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/schedule_parameters.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->vector_points_file_stream.is_open())
      {
        this->estimator->vector_points_file_stream.open(smc_iteration_directory + "/vector_points.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->vector_points_file_stream.is_open())
      {
        for (auto i = this->all_particles[deque_index].particles.begin();
             i!=this->all_particles[deque_index].particles.end();
             ++i)
        {
          (*i)->write_vector_points(this->estimator->vector_variables,
                                    this->estimator->vector_points_file_stream,
                                    NULL);
        }
        //vector_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/vector_points.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->any_points_file_stream.is_open())
      {
        this->estimator->any_points_file_stream.open(smc_iteration_directory + "/any_points.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->any_points_file_stream.is_open())
      {
        for (auto i = this->all_particles[deque_index].particles.begin();
             i!=this->all_particles[deque_index].particles.end();
             ++i)
        {
          (*i)->write_any_points(this->estimator->any_variables,
                                 this->estimator->any_points_file_stream);
        }
        //any_points_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/any_points.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->ancestor_index_file_stream.is_open())
      {
        this->estimator->ancestor_index_file_stream.open(smc_iteration_directory + "/ancestor_index.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->ancestor_index_file_stream.is_open())
      {
        if (deque_index>0)
        {
          for (auto i = this->all_particles[deque_index-1].ancestor_variables.begin();
               i!=this->all_particles[deque_index-1].ancestor_variables.end();
               ++i)
          {
            this->estimator->ancestor_index_file_stream << *i << std::endl;
          }
        }
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/ancestor_index.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->normalised_weights_file_stream.is_open())
      {
        this->estimator->normalised_weights_file_stream.open(smc_iteration_directory + "/normalised_log_weights.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->normalised_weights_file_stream.is_open())
      {
        this->estimator->normalised_weights_file_stream << this->all_particles[deque_index].normalised_log_weights;
        //normalised_weights_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/normalised_log_weights.txt" + " cannot be opened.");
      }
      
      if(!this->estimator->unnormalised_weights_file_stream.is_open())
      {
        this->estimator->unnormalised_weights_file_stream.open(smc_iteration_directory + "/unnormalised_log_weights.txt",std::ios::out | std::ios::app);
      }
      if(this->estimator->unnormalised_weights_file_stream.is_open())
      {
        this->estimator->unnormalised_weights_file_stream << this->all_particles[deque_index].unnormalised_log_weights;
        //unnormalised_weights_file_stream.close();
      }
      else
      {
        Rcpp::stop("File " + smc_iteration_directory + "/unnormalised_log_weights.txt" + " cannot be opened.");
      }
      
      // for any other info in Particles, write to a different file
      
      //for (auto i = this->all_particles[deque_index].particles.begin();
      //     i!=this->all_particles[deque_index].particles.end();
      //     ++i)
      for (size_t i = 0;
           i<this->all_particles[deque_index].particles.size();
           ++i)
      {
        this->all_particles[deque_index].particles[i]->write_factors(smc_iteration_directory,
                                                                     std::to_string(i));
      }
      
      this->close_ofstreams(deque_index);
    }
  }
  
  this->iteration_written_to_file = this->smc_iteration;
  
}

void SMCOutput::close_ofstreams()
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->time_file_stream.close();
  this->estimator->vector_variables_file_stream.close();
  this->estimator->vector_variable_sizes_file_stream.close();
  
  this->estimator->incremental_log_likelihood_file_stream.close();
  this->estimator->output_lengths_file_stream.close();
  this->estimator->resampled_file_stream.close();
  this->estimator->ess_file_stream.close();
  this->estimator->schedule_parameters_file_stream.close();
  this->estimator->vector_points_file_stream.close();
  this->estimator->any_points_file_stream.close(); // should be one for each member of Parameters
  this->estimator->ancestor_index_file_stream.close();
  this->estimator->normalised_weights_file_stream.close();
  this->estimator->unnormalised_weights_file_stream.close();
  
  for (auto i = this->all_particles.begin();
       i!=this->all_particles.end();
       ++i)
  {
    i->close_ofstreams();
  }
}

void SMCOutput::close_ofstreams(size_t deque_index)
{
  this->estimator->log_likelihood_file_stream.close();
  this->estimator->time_file_stream.close();
  this->estimator->vector_variables_file_stream.close();
  this->estimator->vector_variable_sizes_file_stream.close();
  
  this->estimator->incremental_log_likelihood_file_stream.close();
  this->estimator->output_lengths_file_stream.close();
  this->estimator->resampled_file_stream.close();
  this->estimator->ess_file_stream.close();
  this->estimator->schedule_parameters_file_stream.close();
  this->estimator->vector_points_file_stream.close();
  this->estimator->any_points_file_stream.close(); // should be one for each member of Parameters
  this->estimator->ancestor_index_file_stream.close();
  this->estimator->normalised_weights_file_stream.close();
  this->estimator->unnormalised_weights_file_stream.close();
  
  this->all_particles[deque_index].close_ofstreams();
}

void SMCOutput::increment_smc_iteration()
{
  this->smc_iteration = this->smc_iteration + 1;
  //std::cout << this->smc_iteration << std::endl;
}

void SMCOutput::decrement_smc_iteration()
{
  this->smc_iteration = this->smc_iteration - 1;
  //std::cout << this->smc_iteration << std::endl;
}

void SMCOutput::print(std::ostream &os) const
{
  os << "all_particles" << std::endl << "(" << std::endl;
  std::deque<Particles>::const_iterator it;
  for (it=this->all_particles.begin();it!=this->all_particles.end();++it)
  {
    if (it==this->all_particles.begin())
      os << *it;
    else
      os << std::endl << "," << std::endl << *it;
  }
  os << std::endl << ")" << std::endl;

  os << "all_proposed" << std::endl << "(" << std::endl;
  for (it=this->all_proposed.begin();it!=this->all_proposed.end();++it)
  {
    if (it==this->all_proposed.begin())
      os << *it;
    else
      os << std::endl << "," << std::endl << *it;
  }
  os << std::endl << ")" << std::endl;

  /*
  os << "unnormalised_log_weights" << std::endl << "(" << std::endl;
  std::deque<arma::colvec>::const_iterator itd;
  for (itd=this->unnormalised_log_weights.begin();itd!=this->unnormalised_log_weights.end();++itd)
  {
    if (itd==this->unnormalised_log_weights.begin())
      os << *itd;
    else
      os << std::endl << "," << std::endl << *itd;
  }
  os << std::endl << ")" << std::endl;

  os << "normalised_log_weights" << std::endl << "(" << std::endl;
  for (itd=this->normalised_log_weights.begin();itd!=this->normalised_log_weights.end();++itd)
  {
    if (itd==this->normalised_log_weights.begin())
      os << *itd;
    else
      os << std::endl << "," << std::endl << *itd;
  }
  os << std::endl << ")" << std::endl;

  os << "log_normalising_constant_ratios" << std::endl << "(" << std::endl;
  std::vector<double>::const_iterator i;
  for (i=this->log_normalising_constant_ratios.begin();i!=this->log_normalising_constant_ratios.end();++i)
  {
    if (i==this->log_normalising_constant_ratios.begin())
      os << *i;
    else
      os << std::endl << "," << std::endl << *i;
  }
  os << std::endl << ")" << std::endl;
  */
}
