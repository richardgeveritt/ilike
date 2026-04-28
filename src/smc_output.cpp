#include <algorithm>
#include <numeric>
#include <sstream>
#include "smc_output.h"
#include "smc.h"
#include "importance_sampler.h"
#include "smc_mcmc_move.h"
#include "smc_marginal.h"
#include "smc_generic.h"
#include "filesystem.h"
#include "move_output.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{
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
  if (lag_in<2)
  {
    this->lag = 2;
  }
  else
  {
    this->lag = lag_in;
  }
  this->output_lag = lag_in;
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

SMCOutput* SMCOutput::smc_duplicate() const
{
  return( new SMCOutput(*this));
}

void SMCOutput::make_copy(const SMCOutput &another)
{
  this->all_particles = another.all_particles;
  this->all_proposed = another.all_proposed;
  //this->log_normalising_constant_ratios = another.log_normalising_constant_ratios;
  this->lag = another.lag;
  this->output_lag = another.output_lag;
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
  
  if (this->results_name!="")
    this->write(results_name);
}

void SMCOutput::set_time_and_reset_start()
{
  this->set_time();
  
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

void SMCOutput::set_llhd(double llhd_in)
{
  size_t num_to_pop_front = std::max<int>(0,this->llhds.size()-this->lag+1);
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    this->llhds.pop_front();
  }
  this->llhds.push_back(llhd_in);
}

void SMCOutput::forget_you_were_already_written_to_file()
{
  this->iteration_written_to_file = -1;
}

void SMCOutput::terminate()
{
  while (this->all_particles.size()>this->output_lag)
  {
    this->all_particles.pop_back();
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

void SMCOutput::write_to_file(const std::string &dir_name,
                              const std::string &index)
{
  std::string directory_name = dir_name + "_smc";
  
  if (!directory_exists(directory_name))
  {
    make_directory(directory_name);
  }
  
  // Open or create the HDF5 file once per algorithm run
  if (!this->estimator->h5_file)
  {
    this->estimator->h5_file_path = directory_name + "/output.h5";
    this->estimator->h5_file = h5_open_or_create(this->estimator->h5_file_path);
    
    // Store variable names/sizes as root attributes (constant across all iterations)
    auto root = this->estimator->h5_file->getGroup("/");
    h5_set_str_attr(root, "variable_names", this->estimator->vector_variables);
    std::vector<size_t> vsizes(this->estimator->vector_variable_sizes.begin(),
                               this->estimator->vector_variable_sizes.end());
    h5_set_sizet_attr(root, "variable_sizes", vsizes);
  }
  
  HighFive::File &hf = *this->estimator->h5_file;
  
  // for each iteration left to write
  for (size_t iteration = this->iteration_written_to_file+1;
       iteration<this->smc_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->smc_iteration-iteration;
    
    if (this->all_particles.size() > distance_from_end)
    {
      size_t deque_index = this->all_particles.size()-1-distance_from_end;
      
      // ---- top-level extendable datasets ----
      h5_append_double(hf, "log_likelihood", this->llhds[deque_index]);
      h5_append_double(hf, "time",           this->times[deque_index]);
      
      // output_lengths: one row per iteration (one value per particle)
      arma::rowvec ol = this->all_particles[deque_index].get_output_lengths();
      std::vector<double> ol_vec(ol.begin(), ol.end());
      h5_append_row(hf, "output_lengths", ol_vec);
      
      // ---- per-iteration group ----
      std::string iter_grp_path = "iteration/" + std::to_string(iteration + 1);
      HighFive::Group iter_grp = h5_ensure_group(hf, iter_grp_path);
      
      h5_write_scalar_double(iter_grp, "incremental_log_likelihood",
                              this->all_particles[deque_index].log_normalising_constant_ratio);
      h5_write_scalar_int(iter_grp, "resampled",
                           static_cast<int>(this->all_particles[deque_index].resampled_flag));
      h5_write_scalar_double(iter_grp, "ess",
                              this->all_particles[deque_index].ess);
      
      // schedule_parameters as a string
      {
        std::ostringstream oss;
        oss << this->all_particles[deque_index].schedule_parameters;
        h5_write_string(iter_grp, "schedule_parameters", oss.str());
      }
      
      // vector_points: collect all particle rows into one matrix
      {
        std::vector<arma::mat> mats;
        mats.reserve(this->all_particles[deque_index].particles.size());
        size_t total_rows = 0;
        for (auto& p : this->all_particles[deque_index].particles)
        {
          arma::mat m = p->get_matrix_of_vector_points(
                          this->estimator->vector_variables, nullptr);
          total_rows += m.n_rows;
          mats.push_back(std::move(m));
        }
        if (total_rows > 0 && !mats.empty())
        {
          size_t ncols = mats[0].n_cols;
          arma::mat combined(total_rows, ncols);
          size_t row_offset = 0;
          for (auto& m : mats)
          {
            if (m.n_rows > 0)
            {
              combined.rows(row_offset, row_offset + m.n_rows - 1) = m;
              row_offset += m.n_rows;
            }
          }
          h5_write_mat(iter_grp, "vector_points", combined);
        }
      }
      
      // weights
      h5_write_vec(iter_grp, "normalised_log_weights",
                   this->all_particles[deque_index].normalised_log_weights);
      h5_write_vec(iter_grp, "unnormalised_log_weights",
                   this->all_particles[deque_index].unnormalised_log_weights);
      
      // ancestor_index from previous iteration's ancestor_variables
      if (deque_index > 0)
      {
        const auto& av = this->all_particles[deque_index-1].ancestor_variables;
        std::vector<size_t> ai(av.begin(), av.end());
        h5_write_intvec(iter_grp, "ancestor_index", ai);
      }
      else
      {
        h5_write_intvec(iter_grp, "ancestor_index", {});
      }
      
      // nested factor outputs: each creates its own output.h5 in a subdirectory
      std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(iteration+1);
      if (!directory_exists(smc_iteration_directory))
      {
        make_directory(smc_iteration_directory);
      }
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
  // HDF5 file stays open across iterations; reset at algorithm end.
  this->estimator->h5_file.reset();
  
  for (auto i = this->all_particles.begin();
       i!=this->all_particles.end();
       ++i)
  {
    i->close_ofstreams();
  }
}

void SMCOutput::close_ofstreams(size_t deque_index)
{
  // Keep HDF5 file open; only close particles' nested outputs.
  this->all_particles[deque_index].close_ofstreams();
}

void SMCOutput::increment_smc_iteration()
{
  this->smc_iteration = this->smc_iteration + 1;
}

void SMCOutput::decrement_smc_iteration()
{
  this->smc_iteration = this->smc_iteration - 1;
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
}
