#include "transform.h"
#include "ensemble_kalman_output.h"
#include "ensemble_factor_variables.h"
#include "ensemble_kalman.h"
#include "ensemble_kalman_worker.h"
#include "move_output.h"
#include "filesystem.h"
#include "utils.h"
#include "ilike_hdf5_utils.h"
#include <sstream>

namespace ilike
{
EnsembleKalmanOutput::EnsembleKalmanOutput()
:LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
  this->estimator.reset();
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
  this->estimator.reset(estimator_in->ensemble_kalman_duplicate());
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
  this->estimator.reset();
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
  
  this->estimator.reset();
  
  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* EnsembleKalmanOutput::duplicate() const
{
  return( new EnsembleKalmanOutput(*this));
}

void EnsembleKalmanOutput::make_copy(const EnsembleKalmanOutput &another)
{
  this->followed_steps = another.followed_steps;
  this->committed_log_likelihood = another.committed_log_likelihood;
  this->caught_up = another.caught_up;
  if (another.estimator!=NULL)
    this->estimator.reset(another.estimator->ensemble_kalman_duplicate());
  else
    this->estimator.reset();
  this->enk_iteration = another.enk_iteration;
  this->iteration_written_to_file = another.iteration_written_to_file;
  this->transform = another.transform;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
  this->all_ensembles = another.all_ensembles;
  // (C265) the copied ensembles' members point at this output's own copy of
  // the filter, not at the one copied from (deleted with its output: a
  // following filter, which goes on moving a copied output's ensemble, then
  // read freed memory)
  if (this->estimator != NULL)
  {
    EnsembleFactors* factors = this->estimator->ensemble_factors.get();
    for (Ensemble &ensemble : this->all_ensembles)
    {
      ensemble.set_ensemble_factors(factors);
      for (auto &member : ensemble.members)
        if ((member != NULL) && (member->back().ensemble_factor_variables != NULL))
          member->back().ensemble_factor_variables->repoint_ensemble_factors(factors);
      for (auto &member : ensemble.predicted_members)
        if ((member != NULL) && (member->back().ensemble_factor_variables != NULL))
          member->back().ensemble_factor_variables->repoint_ensemble_factors(factors);
    }
  }
  this->lag = another.lag;
  this->output_lag = another.output_lag;
  this->start_time = another.start_time;
  this->times = another.times;
  this->llhds = another.llhds;
  this->skip_to_end_of_sequence = another.skip_to_end_of_sequence;
  this->results_name = another.results_name;
}

Ensemble* EnsembleKalmanOutput::add_ensemble(EnsembleFactors* ensemble_factors)
{
  size_t num_to_pop_front = std::max<int>(0,this->all_ensembles.size()-this->lag+1);
  // (C264) the dropped ensemble's members are kept for the new one's moves
  std::vector<std::unique_ptr<MoveOutput>> recycled;
  for (size_t i=0; i<num_to_pop_front; ++i)
  {
    if (this->estimator->reuses_member_objects() && recycled.empty())
      recycled.swap(this->all_ensembles.front().members);
    this->all_ensembles.pop_front();
  }
  this->all_ensembles.push_back(Ensemble(ensemble_factors));
  this->all_ensembles.back().spare.swap(recycled);
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

double EnsembleKalmanOutput::calculate_mc_inversion_latest_log_normalising_constant_ratio(const Index* index,
                                                                                          double inverse_incremental_temperature)
{
  this->estimator->the_worker->weight(&this->back(),
                                      index,
                                      inverse_incremental_temperature);
  return this->all_ensembles.back().calculate_mc_inversion_log_normalising_constant(this->estimator->the_worker->get_unnormalised_log_incremental_weights(),inverse_incremental_temperature);

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

// (C270) after a parallel chunk group: this output's own copy of the filter
// was made from the group's copy, so it would draw from the group's
// generator (gone with the group); it draws from the main filter's instead
bool EnsembleKalmanOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  if ((estimator_in == NULL) || (this->estimator == NULL))
    return false;
  this->estimator->redirect_random_numbers(estimator_in->get_rng_pointer(), estimator_in->get_seed_pointer());
  return true;
}

// (dev/plan.md 3.21 (c2''')) the caller's move stage
void EnsembleKalmanOutput::advance_at_caller_target(const Parameters &parameters)
{
  this->estimator->advance_at_caller_target_enk(this, parameters);
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

void EnsembleKalmanOutput::subsample_evaluate_smcfixed_part(const Parameters &conditioned_on_parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->ensemble_kalman_subsample_evaluate(this, conditioned_on_parameters);
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
  return this->estimator.get();
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

// the output file (made with its attributes the first time)
HighFive::File& EnsembleKalmanOutput::output_file(const std::string &dir_name)
{
  std::string directory_name = dir_name + "_enk";

  if (!directory_exists(directory_name))
    make_directory(directory_name);

  if (!this->estimator->h5_file)
  {
    this->estimator->h5_file_path = directory_name + "/output.h5";
    this->estimator->h5_file = h5_open_or_create(this->estimator->h5_file_path);

    auto root = this->estimator->h5_file->getGroup("/");
    h5_set_str_attr(root, "variable_names", this->estimator->vector_variables);
    std::vector<size_t> vsizes(this->estimator->vector_variable_sizes.begin(),
                               this->estimator->vector_variable_sizes.end());
    h5_set_sizet_attr(root, "variable_sizes", vsizes);
  }
  return *this->estimator->h5_file;
}

// (C264) with nothing stored, write_to_file makes the file and writes no
// iteration: the same here
void EnsembleKalmanOutput::write_compact_values(const std::string &dir_name,
                                                const std::string &index,
                                                const double *values)
{
  this->output_file(dir_name);
}

void EnsembleKalmanOutput::write_to_file(const std::string &dir_name,
                                         const std::string &index)
{
  HighFive::File &hf = this->output_file(dir_name);
  std::string directory_name = dir_name + "_enk";

  for (size_t iteration = this->iteration_written_to_file+1;
       iteration < this->enk_iteration+1;
       ++iteration)
  {
    size_t distance_from_end = this->enk_iteration - iteration;

    if (this->all_ensembles.size() > distance_from_end)
    {
      size_t deque_index = this->all_ensembles.size()-1-distance_from_end;

      h5_append_double(hf, "log_likelihood", this->llhds[deque_index]);
      h5_append_double(hf, "time",           this->times[deque_index]);

      {
        arma::rowvec ol = this->all_ensembles[deque_index].get_output_lengths();
        std::vector<double> ol_vec(ol.begin(), ol.end());
        h5_append_row(hf, "output_lengths", ol_vec);
      }

      std::string iter_grp_path = "iteration/" + std::to_string(iteration+1);
      HighFive::Group iter_grp = h5_ensure_group(hf, iter_grp_path);

      h5_write_scalar_double(iter_grp, "incremental_log_likelihood",
                             this->all_ensembles[deque_index].log_normalising_constant_ratio);
      h5_write_scalar_double(iter_grp, "ess",
                             this->all_ensembles[deque_index].ess);
      {
        std::ostringstream oss;
        if (this->estimator->reciprocal_schedule_scale == 0.0)
          oss << this->all_ensembles[deque_index].schedule_parameters;
        else
          oss << this->estimator->reciprocal_schedule_scale *
                 pow(this->all_ensembles[deque_index].schedule_parameters, -0.5);
        h5_write_string(iter_grp, "schedule_parameters", oss.str());
      }

      // vector_points: collect all ensemble members into one matrix (C279: a
      // compact ensemble's rows are that matrix: its members hold only the
      // state variables, without a transform)
      if (this->all_ensembles[deque_index].compact)
      {
        const arma::mat &rows = this->all_ensembles[deque_index].packed_members;
        if ((rows.n_rows > 0) && (this->transform == NULL))
          h5_write_mat(iter_grp, "vector_points", rows);
        else if (rows.n_rows > 0)
        {
          // (C281, MFDS) each member through the inverse transform, as a
          // member object's get_matrix_of_vector_points
          const PackingInstructions &packing = this->estimator->packing_instructions;
          arma::mat points;
          for (arma::uword i = 0; i < rows.n_rows; ++i)
          {
            Parameters member;
            for (size_t j = 0; j < packing.states_names.size(); ++j)
              member[packing.states_names[j]] = arma::colvec(rows.row(i).subvec(packing.states_start_and_end[j].first,
                                                                                 packing.states_start_and_end[j].second).t());
            arma::rowvec point = this->transform->inverse_transform(member).get_rowvec(this->estimator->vector_variables);
            if (i == 0)
              points.set_size(rows.n_rows, point.n_elem);
            points.row(i) = point;
          }
          h5_write_mat(iter_grp, "vector_points", points);
        }
      }
      else
      {
        std::vector<arma::mat> mats;
        mats.reserve(this->all_ensembles[deque_index].members.size());
        size_t total_rows = 0;
        for (auto &m : this->all_ensembles[deque_index].members)
        {
          arma::mat row_mat = m->get_matrix_of_vector_points(
                                this->estimator->vector_variables, this->transform);
          total_rows += row_mat.n_rows;
          mats.push_back(std::move(row_mat));
        }
        if (total_rows > 0 && !mats.empty())
        {
          size_t ncols = mats[0].n_cols;
          arma::mat combined(total_rows, ncols);
          size_t row_offset = 0;
          for (auto &m : mats)
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

      std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(iteration+1);
      if (!directory_exists(smc_iteration_directory))
        make_directory(smc_iteration_directory);

      for (size_t i = 0; i < this->all_ensembles[deque_index].members.size(); ++i)
        this->all_ensembles[deque_index].members[i]->write_factors(smc_iteration_directory,
                                                                   std::to_string(i));

      for (size_t i = 0; i < this->all_ensembles[deque_index].members.size(); ++i)
        this->all_ensembles[deque_index].members[i]->write_ensemble_factors(smc_iteration_directory,
                                                                            std::to_string(i));

      this->close_ofstreams(deque_index);
    }
  }

  this->iteration_written_to_file = this->enk_iteration;
}

void EnsembleKalmanOutput::close_ofstreams()
{
  this->estimator->h5_file.reset();

  for (auto i = this->all_ensembles.begin();
       i != this->all_ensembles.end();
       ++i)
  {
    i->close_ofstreams();
  }
}

void EnsembleKalmanOutput::close_ofstreams(size_t deque_index)
{
  this->all_ensembles[deque_index].close_ofstreams();
}

// Adaptive summary (dev/plan.md 2.6a): only when smcfixed, where the adaptive
// part leaves log_likelihood as the smcfixed part set it.
int EnsembleKalmanOutput::adaptive_summary_size() const
{
  return this->estimator->smcfixed_flag ? 1 : -1;
}

void EnsembleKalmanOutput::adaptive_summary(double *values) const
{
  values[0] = this->log_likelihood;
}

double EnsembleKalmanOutput::smcadaptive_from_summary(const double *values,
                                     const Parameters &parameters) const
{
  return values[0];
}

}
