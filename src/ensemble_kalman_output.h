#ifndef ENSEMBLEKALMANOUTPUT_H
#define ENSEMBLEKALMANOUTPUT_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <deque>
#include <chrono>

#include "likelihood_estimator_output.h"
#include "ilike_hdf5_utils.h"
#include "ensemble.h"
#include "ilike_header.h"

namespace ilike
{

class EnsembleKalman;
class EnsembleKalmanFilter;
class EnsembleSequencer;
class EnsembleKalmanInversion;
class EnsembleKalmanMFDS;

class EnsembleKalmanOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  // unsure
  // copied from KF
  
  EnsembleKalmanOutput();
  EnsembleKalmanOutput(EnsembleKalman* estimator_in,
                       size_t lag_in,
                       std::shared_ptr<Transform> transform_in,
                       const std::string &results_name_in);
  virtual ~EnsembleKalmanOutput();
  
  EnsembleKalmanOutput(const EnsembleKalmanOutput &another);
  void operator=(const EnsembleKalmanOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  
  Ensemble* add_ensemble(EnsembleFactors* ensemble_factors);
  void initialise_next_step();
  
  Ensemble back() const;
  Ensemble& back();
  
  double calculate_latest_log_normalising_constant_ratio();
  double calculate_mc_inversion_latest_log_normalising_constant_ratio(const Index* index,
                                                                      double inverse_incremental_temperature);
  double calculate_inversion_latest_log_normalising_constant_ratio(double inverse_incremental_temperature);
  double calculate_unbiased_inversion_latest_log_normalising_constant_ratio(double inverse_incremental_temperature);
  double calculate_path1_inversion_latest_log_normalising_constant_ratio(const std::vector<double> &previous_log_measurement_likelihood_means,
                                                                         double inverse_incremental_temperature,
                                                                         double temperature,
                                                                         double multiplier);
  double calculate_path2_inversion_latest_log_normalising_constant_ratio(const std::vector<double> &previous_log_measurement_likelihood_means,
                                                                         const std::vector<double> &previous_log_measurement_likelihood_variances,
                                                                         double inverse_incremental_temperature);
  
  void calculate_path1_inversion_initial_quantities(std::vector<double> &initial_log_measurement_likelihood_means,
                                                    double temperature,
                                                    double multiplier);
  void calculate_path2_inversion_initial_quantities(std::vector<double> &initial_log_measurement_likelihood_means,
                                                    std::vector<double> &initial_log_measurement_likelihood_variances);
  
  void calculate_kalman_gains(double inverse_incremental_temperature);
  
  void simulate();
  
  void simulate(const Parameters &parameters);

  void evaluate_smcfixed_part(const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters);
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;

  void subsample_simulate();
  void subsample_simulate(const Parameters &parameters);

  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  size_t number_of_ensemble_kalman_iterations() const;
  
  void increment_enk_iteration();
  
  void set_time();
  
  void set_llhd(double llhd_in);

  void forget_you_were_already_written_to_file();
  
  void terminate();
  
  void skip_to_end_of_sequence_if_points_are_gaussian(double significance_level);
  
  // (C264) a compact record when nothing is stored (iterations_to_store =
  // 0, the default inside another algorithm): the log-likelihood, and the
  // file write_to_file would make (its attributes only); an MCMC chain then
  // stores each iteration as numbers instead of a copy of the particle
  size_t compact_size() const { return (this->output_lag == 0) ? 1 : 0; }
  void compact_values(double *values) const { values[0] = this->log_likelihood; }
  void write_compact_values(const std::string &directory_name,
                            const std::string &index,
                            const double *values);
  HighFive::File& output_file(const std::string &dir_name);
  // (dev/plan.md 3.21 (c2''')) following the caller's data tempering: the
  // observations committed, their log-likelihood, and whether this output has
  // been run through the caller's committed observations yet
  size_t followed_steps = 0;
  double committed_log_likelihood = 0.0;
  bool caught_up = false;
  void advance_at_caller_target(const Parameters &parameters);
  void close_ofstreams();
  
protected:
  
  friend EnsembleKalmanFilter;
  friend EnsembleKalmanInversion;
  friend EnsembleSequencer;
  friend EnsembleKalmanMFDS;
  
  // Stored in Factors.
  // Moving to be stored here (May 2024).
  std::unique_ptr<EnsembleKalman> estimator;
  
  size_t enk_iteration;
  
  int iteration_written_to_file;
  
  std::string results_name;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;
  
  bool skip_to_end_of_sequence;
  
  std::shared_ptr<Transform> transform;
  
  void close_ofstreams(size_t deque_index);
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  void make_copy(const EnsembleKalmanOutput &another);
  
  std::deque<Ensemble> all_ensembles;
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
  size_t lag;
  size_t output_lag;
  
  std::chrono::high_resolution_clock::time_point start_time;
  
};
}

#endif
