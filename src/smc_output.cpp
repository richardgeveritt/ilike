#include <algorithm>
#include <numeric>
#include <sstream>
#include "smc_output.h"
#include "factor_variables.h"
#include "utils.h"
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
  this->pending = another.pending;
  this->all_particles = another.all_particles;
  this->all_proposed = another.all_proposed;
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
  this->compact_generations = another.compact_generations;
  this->compacted = another.compacted;
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
  // (dev/plan.md 3.21) deferred to evaluation when the estimator allows it
  if (this->estimator->defers_simulation())
  {
    this->pending = true;
    return;
  }
  this->estimator->simulate_smc(this, parameters);
}

void SMCOutput::run_if_pending(const Parameters &parameters)
{
  if (!this->pending)
    return;
  this->pending = false;
  this->estimator->simulate_smc(this, parameters);
}

double SMCOutput::evaluate_with_threshold(const Parameters &parameters, double threshold, bool &below)
{
  below = false;
  if (!this->pending)
    return this->evaluate(parameters);
  this->pending = false;
  if (this->estimator->simulate_smc_with_threshold(this, parameters, threshold))
  {
    below = true;
    return -arma::datum::inf;
  }
  return this->evaluate(parameters);
}

void SMCOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  this->run_if_pending(parameters);
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smc(this, parameters);
  }
  else
  {
    this->estimator->evaluate_smcfixed_part_smc(this, parameters);
  }
}

void SMCOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  this->run_if_pending(parameters);
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
  // the objects of the dropped generation are kept for the new one to reuse
  // when the algorithm can (SMC::reuses_particle_objects)
  std::vector<std::unique_ptr<MoveOutput>> recycled;
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    if (this->estimator->reuses_particle_objects() && recycled.empty())
      recycled.swap(this->all_particles.front().particles);
    this->all_particles.pop_front();
  }
  this->all_particles.push_back(Particles());
  this->all_particles.back().spare.swap(recycled);
  this->all_particles.back().reserve(this->estimator->number_of_particles);
  return &this->all_particles.back();
}

Particles* SMCOutput::add_particles(Particles* most_recent_particles)
{
  size_t num_to_pop_front = std::max<int>(0,this->all_particles.size()-this->lag+1);
  // the objects of the dropped generation are kept for the new one to reuse
  // when the algorithm can (SMC::reuses_particle_objects)
  std::vector<std::unique_ptr<MoveOutput>> recycled;
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    if (this->estimator->reuses_particle_objects() && recycled.empty())
      recycled.swap(this->all_particles.front().particles);
    this->all_particles.pop_front();
  }
  this->all_particles.push_back(Particles());
  this->all_particles.back().spare.swap(recycled);
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

void SMCOutput::count_divergences(size_t &divergent, size_t &iterations) const
{
  divergent = 0;
  iterations = 0;
  if (this->all_particles.empty())
    return;
  for (auto &p : this->all_particles.back().particles)
    if (p)
      p->count_divergences(divergent, iterations);
}

void SMCOutput::terminate()
{
  // Keep the output_lag most recent generations. (These used to be removed
  // from the back, i.e. the newest ones: with one generation to store, the
  // final generation was dropped and the one before was written as if it
  // were the final one, B24.)
  while (this->all_particles.size()>this->output_lag)
  {
    this->all_particles.pop_front();
  }
  
  while (this->llhds.size()>this->output_lag)
  {
    this->llhds.pop_front();
  }
  
  while (this->times.size()>this->output_lag)
  {
    this->times.pop_front();
  }

  if (this->estimator->compact_stored_output)
    this->compact();
}

void SMCOutput::compact()
{
  this->compact_generations.clear();
  size_t number_stored = this->all_particles.size();
  size_t ncols = 0;
  for (auto v : this->estimator->vector_variable_sizes)
    ncols += v;
  for (size_t g = 0; g < number_stored; ++g)
  {
    Particles &particles = this->all_particles[g];
    CompactGeneration generation;
    generation.iteration = this->smc_iteration - (number_stored - 1 - g);
    arma::rowvec ol = particles.get_output_lengths();
    generation.output_lengths.assign(ol.begin(), ol.end());
    size_t total = size_t(arma::accu(ol));
    generation.points_transposed.set_size(ncols, total);
    size_t column = 0;
    if (particles.is_flat() && (total > 0))
      particles.flat_points.fill_transposed(this->estimator->vector_variables,
                                            generation.points_transposed, 0);
    for (auto &p : particles.particles)
      column += p->fill_vector_points_transposed(this->estimator->vector_variables,
                                                 generation.points_transposed, column);
    generation.unnormalised_log_weights = particles.unnormalised_log_weights;
    generation.normalised_log_weights = particles.normalised_log_weights;
    if (g > 0)
      generation.ancestors = this->all_particles[g-1].ancestor_variables;
    generation.incremental_log_likelihood = particles.log_normalising_constant_ratio;
    generation.ess = particles.ess;
    generation.resampled = static_cast<int>(particles.resampled_flag);
    std::ostringstream oss;
    oss << particles.schedule_parameters;
    generation.schedule_parameters = oss.str();
    this->compact_generations.push_back(std::move(generation));
  }
  // the rich particles are no longer needed
  this->all_particles.clear();
  this->compacted = true;
}

void SMCOutput::write_compact_nested(const std::string &dir_name)
{
  // One file per outer iteration (directory dir_name + "_smc"), shared by the
  // outer points written one after another: each stored inner generation's
  // points, weights, ancestors and output lengths are appended, so they are
  // concatenated over the outer points; log_likelihood, time and the
  // per-generation scalars have one value per outer point.
  std::string directory_name = dir_name + "_smc";
  if (!this->estimator->h5_file)
  {
    if (!directory_exists(directory_name))
      make_directory(directory_name);
    this->estimator->h5_file_path = directory_name + "/output.h5";
    this->estimator->h5_file = h5_open_or_create(this->estimator->h5_file_path);
    auto root = this->estimator->h5_file->getGroup("/");
    if (!root.hasAttribute("variable_names"))
    {
      h5_set_str_attr(root, "variable_names", this->estimator->vector_variables);
      std::vector<size_t> vsizes(this->estimator->vector_variable_sizes.begin(),
                                 this->estimator->vector_variable_sizes.end());
      h5_set_sizet_attr(root, "variable_sizes", vsizes);
      h5_set_sizet_attr(root, "nested_external_points", {1});
    }
  }
  HighFive::File &hf = *this->estimator->h5_file;
  h5_append_double(hf, "log_likelihood", this->log_likelihood);
  h5_append_double(hf, "time", this->times.empty() ? 0.0 : this->times.back());
  for (auto &generation : this->compact_generations)
  {
    std::string group = "iteration/" + std::to_string(generation.iteration + 1);
    HighFive::Group grp = h5_ensure_group(hf, group);
    h5_append_rows(hf, group + "/vector_points", generation.points_transposed.memptr(),
                   generation.points_transposed.n_cols, generation.points_transposed.n_rows);
    h5_append_doubles(hf, group + "/unnormalised_log_weights",
                      std::vector<double>(generation.unnormalised_log_weights.begin(), generation.unnormalised_log_weights.end()));
    h5_append_doubles(hf, group + "/normalised_log_weights",
                      std::vector<double>(generation.normalised_log_weights.begin(), generation.normalised_log_weights.end()));
    if (!generation.ancestors.empty())
      h5_append_doubles(hf, group + "/ancestor_index",
                        std::vector<double>(generation.ancestors.begin(), generation.ancestors.end()));
    h5_append_doubles(hf, group + "/output_lengths", generation.output_lengths);
    h5_append_double(hf, group + "/incremental_log_likelihood", generation.incremental_log_likelihood);
    h5_append_double(hf, group + "/ess", generation.ess);
    h5_append_int(hf, group + "/resampled", generation.resampled);
    if (!grp.exist("schedule_parameters"))
      h5_write_string(grp, "schedule_parameters", generation.schedule_parameters);
  }
}

void SMCOutput::write_to_file(const std::string &dir_name,
                              const std::string &index)
{
  // A nested run kept in compact form (SMC::compact_stored_output).
  if (this->compacted)
  {
    if (!this->compact_generations.empty())
      this->write_compact_nested(dir_name);
    return;
  }

  // An output set to store no iterations (e.g. the importance sampler inside
  // an ABC likelihood, built without store_output) writes nothing. (Such an
  // output may still hold particles: an ABC likelihood whose tolerance is set
  // by an outer SMC keeps its simulations to reweight them.)
  if (this->output_lag == 0)
    return;

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

      // all_particles, llhds and times hold one entry per stored SMC
      // iteration and must stay aligned.
      if ((deque_index >= this->llhds.size()) || (deque_index >= this->times.size()))
      {
        Rcpp::stop("SMCOutput::write_to_file: output for SMC iteration " +
                   std::to_string(iteration+1) + " is incomplete (" +
                   std::to_string(this->all_particles.size()) + " particle sets, " +
                   std::to_string(this->llhds.size()) + " log-likelihoods, " +
                   std::to_string(this->times.size()) + " times).");
      }

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
      
      // vector_points: all particle rows in one matrix. Its size is known
      // (rows: the summed output lengths; columns: the summed variable sizes),
      // so each particle's rows are copied straight into it (the matrices
      // used to be collected in a vector first and then copied again).
      {
        size_t total_rows = size_t(arma::accu(ol));
        size_t ncols = 0;
        for (auto s : this->estimator->vector_variable_sizes)
          ncols += s;
        if (total_rows > 0 && ncols > 0)
        {
          // one column per point; written as a (points x columns) dataset
          arma::mat combined_transposed(ncols, total_rows);
          size_t column_offset = 0;
          if (this->all_particles[deque_index].is_flat())
          {
            this->all_particles[deque_index].flat_points.fill_transposed(this->estimator->vector_variables,
                                                                         combined_transposed,
                                                                         0);
            column_offset = total_rows;
          }
          for (auto& p : this->all_particles[deque_index].particles)
          {
            size_t remaining = total_rows - column_offset;
            if (size_t(p->length()) > remaining)
              Rcpp::stop("SMCOutput::write_to_file: a particle's points do not match the expected output size.");
            column_offset += p->fill_vector_points_transposed(this->estimator->vector_variables,
                                                              combined_transposed,
                                                              column_offset);
          }
          if (column_offset != total_rows)
            Rcpp::stop("SMCOutput::write_to_file: the particles' points do not match the stored output lengths.");
          h5_write_mat_transposed(iter_grp, "vector_points", combined_transposed);
        }
      }
      
      // adaptation output of adapting MCMC moves (dev/plan.md 3.2): a trace
      // row for each stored point (aligned with vector_points) and the
      // state at the end of each chain; column names as root attributes
      {
        const auto &particles = this->all_particles[deque_index].particles;
        const MoveOutput* first_adapting = NULL;
        for (auto &p : particles)
        {
          if (p->adaptation_trace_columns() > 0)
          {
            first_adapting = p.get();
            break;
          }
        }
        if ( (first_adapting!=NULL) && !this->all_particles[deque_index].is_flat() )
        {
          size_t trace_columns = first_adapting->adaptation_trace_columns();
          size_t final_columns = first_adapting->adaptation_final_columns();
          size_t total_rows = 0;
          for (auto &p : particles)
            total_rows += p->length();
          arma::mat trace_transposed(trace_columns, total_rows, arma::fill::zeros);
          arma::mat final_transposed(final_columns, particles.size(), arma::fill::zeros);
          size_t column_offset = 0;
          for (size_t i=0; i<particles.size(); ++i)
          {
            if (particles[i]->adaptation_trace_columns()!=trace_columns)
              Rcpp::stop("SMCOutput::write_to_file: not every chain has the same adaptation output.");
            size_t rows = particles[i]->fill_adaptation_trace_transposed(trace_transposed, column_offset);
            if (rows!=particles[i]->length())
              Rcpp::stop("SMCOutput::write_to_file: a chain's adaptation trace does not match its stored points.");
            column_offset += rows;
            particles[i]->fill_adaptation_final(final_transposed, i);
          }
          h5_write_mat_transposed(iter_grp, "adaptation_trace", trace_transposed);
          h5_write_mat_transposed(iter_grp, "adaptation_final", final_transposed);
          std::vector<std::string> trace_names, final_names;
          first_adapting->adaptation_names(trace_names, final_names);
          auto root = hf.getGroup("/");
          h5_set_str_attr(root, "adaptation_trace_names", trace_names);
          h5_set_str_attr(root, "adaptation_final_names", final_names);
        }
      }

      // per-iteration diagnostics (nuts; C212): a row for each stored point,
      // aligned with vector_points
      {
        const auto &particles = this->all_particles[deque_index].particles;
        const MoveOutput* first = NULL;
        for (auto &p : particles)
        {
          if (p->diagnostics_columns() > 0)
          {
            first = p.get();
            break;
          }
        }
        if ( (first!=NULL) && !this->all_particles[deque_index].is_flat() )
        {
          size_t columns = first->diagnostics_columns();
          size_t total_rows = 0;
          for (auto &p : particles)
            total_rows += p->length();
          arma::mat transposed(columns, total_rows, arma::fill::zeros);
          size_t column_offset = 0;
          for (size_t i=0; i<particles.size(); ++i)
          {
            if (particles[i]->diagnostics_columns()!=columns)
              Rcpp::stop("SMCOutput::write_to_file: not every chain has the same diagnostics.");
            size_t rows = particles[i]->fill_diagnostics_transposed(transposed, column_offset);
            if (rows!=particles[i]->length())
              Rcpp::stop("SMCOutput::write_to_file: a chain's diagnostics do not match its stored points.");
            column_offset += rows;
          }
          h5_write_mat_transposed(iter_grp, "diagnostics", transposed);
          std::vector<std::string> names;
          first->diagnostics_names(names);
          auto root = hf.getGroup("/");
          h5_set_str_attr(root, "diagnostics_names", names);
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
      // (skipped when the algorithm was asked not to write them)
      if (this->estimator->write_factor_output)
      {
      std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(iteration+1);
      if (!directory_exists(smc_iteration_directory))
      {
        make_directory(smc_iteration_directory);
      }
      Particles &generation = this->all_particles[deque_index];
      if (generation.is_flat())
      {
        // as write_factors would have written for each particle
        if (generation.flat_record_writer!=NULL)
        {
          // each particle's stored points in turn (one for importance
          // sampling, a chain's stored states after MCMC moves), all under
          // the particle's index
          arma::rowvec lengths = generation.get_output_lengths();
          size_t row = 0;
          for (size_t i = 0; i<generation.size(); ++i)
          {
            std::string particle_index = std::to_string(i);
            for (size_t j = 0; j<size_t(lengths[i]); ++j, ++row)
              generation.flat_record_writer->write_compact_values(smc_iteration_directory,
                                                                  particle_index,
                                                                  generation.flat_points.extras_of(row));
          }
        }
        else if (generation.flat_factor_variables.size()==generation.size())
        {
          for (size_t i = 0; i<generation.size(); ++i)
            generation.flat_factor_variables[i]->write_to_file(smc_iteration_directory,
                                                              std::to_string(i));
        }
        else if (generation.size()>0)
        {
          Rcpp::stop("SMCOutput::write_to_file: the factor outputs of this generation were not kept.");
        }
      }
      for (size_t i = 0;
           i<generation.particles.size();
           ++i)
      {
        generation.particles[i]->write_factors(smc_iteration_directory,
                                               std::to_string(i));
      }
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
  
}

// Adaptive summary (dev/plan.md 2.6a): only when smcfixed, where the adaptive
// part leaves log_likelihood as the smcfixed part set it.
int SMCOutput::adaptive_summary_size() const
{
  // smcfixed: the adaptive part leaves log_likelihood as the smcfixed part
  // set it. An importance sampler reweighted at the caller's target (the ABC
  // likelihood under a tolerance sequence): its single flat generation's
  // summaries — the number of values per point s, then each point's s
  // values, its previous target and its previous normalised log weight —
  // from which smcadaptive_from_summary computes the log normalising
  // constant at a trial target exactly as the reweighting would
  // (ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc with
  // weight_from_summaries_if_possible, update_weights and
  // calculate_log_normalising_constant).
  if (this->estimator->smcfixed_flag)
    return 1;
  if (!this->estimator->reweights_at_caller_target() || (this->all_particles.size() != 1))
    return -1;
  const Particles &generation = this->all_particles.back();
  size_t N = generation.size();
  if (!generation.is_flat() || (generation.adaptive_summary_state != 1) ||
      (generation.adaptive_summary_evaluator == NULL) ||
      (generation.adaptive_previous_targets.n_elem != N) ||
      (generation.previous_normalised_log_weights.n_elem != N))
    return -1;
  return int(1 + generation.adaptive_summaries.n_rows*N + 2*N);
}

void SMCOutput::adaptive_summary(double *values) const
{
  if (this->estimator->smcfixed_flag)
  {
    values[0] = this->log_likelihood;
    return;
  }
  const Particles &generation = this->all_particles.back();
  size_t N = generation.size();
  size_t s = generation.adaptive_summaries.n_rows;
  values[0] = double(s);
  std::copy(generation.adaptive_summaries.memptr(), generation.adaptive_summaries.memptr() + s*N, values + 1);
  std::copy(generation.adaptive_previous_targets.memptr(), generation.adaptive_previous_targets.memptr() + N, values + 1 + s*N);
  std::copy(generation.previous_normalised_log_weights.memptr(), generation.previous_normalised_log_weights.memptr() + N, values + 1 + s*N + N);
}

double SMCOutput::smcadaptive_from_summary(const double *values,
                                           const Parameters &parameters) const
{
  if (this->estimator->smcfixed_flag)
    return values[0];
  // this output's own generation supplies the evaluator (factor variables of
  // the inner sampler's factors) and the index; the numbers come from values
  const Particles &generation = this->all_particles.back();
  const FactorVariables* evaluator = generation.adaptive_summary_evaluator;
  const Index* index = generation.adaptive_summary_index;
  size_t N = generation.size();
  size_t s = size_t(values[0]);
  const double *summaries = values + 1;
  const double *previous_targets = values + 1 + s*N;
  const double *previous_weights = values + 1 + s*N + N;
  arma::colvec weights(N);
  std::vector<double> targets(N);
  if (N > 0)
    evaluator->smcadaptive_from_summaries(index, summaries, s, N, parameters, targets.data());
  for (size_t i = 0; i < N; ++i)
  {
    double a = targets[i];
    if (a==-arma::datum::inf)
      weights[i] = -arma::datum::inf;
    else
      weights[i] = a - previous_targets[i];
  }
  // as Particles::update_weights and calculate_log_normalising_constant
  arma::colvec previous(const_cast<double*>(previous_weights), N, false, true);
  arma::colvec unnormalised = previous + weights;
  return log_sum_exp(unnormalised);
}

bool SMCOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  // the stored generations must hold no rich particles (they would point to
  // the copy's factors): compacted, flat, or none
  SMC* smc = dynamic_cast<SMC*>(estimator_in);
  if (smc == NULL)
    return false;
  for (auto g = this->all_particles.begin(); g != this->all_particles.end(); ++g)
  {
    if (!g->repoint_flat_to(smc->factors.get(), smc->get_weighting_index()))
      return false;
  }
  for (auto g = this->all_proposed.begin(); g != this->all_proposed.end(); ++g)
  {
    if (g->size() > 0)
      return false;
  }
  this->estimator = smc;
  return true;
}

}
