#ifndef ENSEMBLEKALMANOUTPUT_H
#define ENSEMBLEKALMANOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>
#include <chrono>

#include "likelihood_estimator_output.h"
#include "ensemble.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file ensemble_kalman_output.h
   * @brief Defines the EnsembleKalman class.
   *
   * An Ensemble Kalman method implementation. Performs ensemble-based data assimilation using Kalman-type updates.
   *
   * @namespace ilike
   * @class EnsembleKalman
   * @brief The ensemble kalman class.
   */


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
  
  /**
   * @brief Performs the ensemblekalmanoutput operation.
   */
  EnsembleKalmanOutput();
  EnsembleKalmanOutput(EnsembleKalman* estimator_in,
                       size_t lag_in,
                       std::shared_ptr<Transform> transform_in,
                       const std::string &results_name_in);
  /**
   * @brief Performs the ~ensemblekalmanoutput operation.
   */
  virtual ~EnsembleKalmanOutput();
  
  /**
   * @brief Performs the ensemblekalmanoutput operation.
   *
   * @param another The EnsembleKalman instance to copy from.
   */
  EnsembleKalmanOutput(const EnsembleKalmanOutput &another);
  /**
   * @brief Assignment operator for EnsembleKalman.
   *
   * @param another The EnsembleKalman instance to copy from.
   */
  void operator=(const EnsembleKalmanOutput &another);
  /**
   * @brief Creates a deep copy of this EnsembleKalman object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* duplicate() const;
  
  //Ensemble* add_ensemble();
  /**
   * @brief Performs the add ensemble operation.
   *
   * @param ensemble_factors The ensemble factors.
   *
   * @return The result.
   */
  Ensemble* add_ensemble(EnsembleFactors* ensemble_factors);
  //void add_proposed_ensemble(const Ensemble &latest_proposals);
  //void initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  /**
   * @brief Performs the initialise next step operation.
   */
  void initialise_next_step();
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Ensemble back() const;
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  Ensemble& back();
  
  /**
   * @brief Performs the calculate latest log normalising constant ratio operation.
   *
   * @return The result.
   */
  double calculate_latest_log_normalising_constant_ratio();
  double calculate_mc_inversion_latest_log_normalising_constant_ratio(const Index* index,
                                                                      double inverse_incremental_temperature);
  /**
   * @brief Performs the calculate inversion latest log normalising constant ratio operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  double calculate_inversion_latest_log_normalising_constant_ratio(double inverse_incremental_temperature);
  /**
   * @brief Performs the calculate unbiased inversion latest log normalising constant ratio operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
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
  
  /**
   * @brief Performs the calculate kalman gains operation.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   */
  void calculate_kalman_gains(double inverse_incremental_temperature);
  
  /**
   * @brief Simulates the required variables.
   */
  void simulate();
  
  /**
   * @brief Simulates the required variables.
   *
   * @param parameters The parameters.
   */
  void simulate(const Parameters &parameters);
  
  //double evaluate(const Parameters &parameters);
  
  /*
   void evaluate_smcfixed_part();
   void evaluate_smcadaptive_part_given_smcfixed();
   
   void subsample_evaluate_smcfixed_part();
   void subsample_evaluate_smcadaptive_part_given_smcfixed();
   */
  
  
  void evaluate_smcfixed_part(const Parameters &conditioned_on_parameters);
  /**
   * @brief Evaluates the smcadaptive part given smcfixed.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters);
  
  
  /**
   * @brief Performs the subsample simulate operation.
   */
  void subsample_simulate();
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param parameters The parameters.
   */
  void subsample_simulate(const Parameters &parameters);
  
  
  /**
   * @brief Performs the subsample evaluate smcfixed part operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  
  /**
   * @brief Returns the likelihood estimator.
   *
   * @return The result.
   */
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  /**
   * @brief Performs the number of ensemble kalman iterations operation.
   *
   * @return The result.
   */
  size_t number_of_ensemble_kalman_iterations() const;
  
  /**
   * @brief Performs the increment enk iteration operation.
   */
  void increment_enk_iteration();
  
  /**
   * @brief Sets the time.
   */
  void set_time();
  
  /**
   * @brief Sets the llhd.
   *
   * @param llhd_in The llhd.
   */
  void set_llhd(double llhd_in);
  
  /*
   void set_current_predicted_statistics(const arma::colvec &latest_mean,
   const arma::mat &latest_covariance);
   void add_predicted_statistics();
   
   void set_current_posterior_statistics(const arma::colvec &latest_mean,
   const arma::mat &latest_covariance);
   void add_posterior_statistics();
   
   
   arma::colvec predicted_mean_back() const;
   arma::colvec posterior_mean_back() const;
   arma::mat predicted_covariance_back() const;
   arma::mat posterior_covariance_back() const;
   
   void set_current_predicted_to_be_current_posterior();
   
   size_t predicted_size() const;
   
   void print(std::ostream &os) const;
   */
  
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   */
  void terminate();
  
  /**
   * @brief Performs the skip to end of sequence if points are gaussian operation.
   *
   * @param significance_level The significance level.
   */
  void skip_to_end_of_sequence_if_points_are_gaussian(double significance_level);
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  friend EnsembleKalmanFilter;
  friend EnsembleKalmanInversion;
  friend EnsembleSequencer;
  friend EnsembleKalmanMFDS;
  
  // Stored in Factors.
  // Moving to be stored here (May 2024).
  /** @brief The estimator. */
  EnsembleKalman* estimator;
  
  /** @brief The enk iteration. */
  size_t enk_iteration;
  
  /** @brief The iteration written to file. */
  int iteration_written_to_file;
  
  /** @brief The results name. */
  std::string results_name;
  
  /** @brief The log likelihood smcfixed part. */
  double log_likelihood_smcfixed_part;
  /** @brief The subsample log likelihood smcfixed part. */
  double subsample_log_likelihood_smcfixed_part;
  
  /** @brief The skip to end of sequence. */
  bool skip_to_end_of_sequence;
  
  /** @brief The transform. */
  std::shared_ptr<Transform> transform;
  
  /**
   * @brief Closes any open file streams.
   *
   * @param deque_index The deque index.
   */
  void close_ofstreams(size_t deque_index);
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  /**
   * @brief Copies the state of another EnsembleKalman into this object.
   *
   * @param another The EnsembleKalman instance to copy from.
   */
  void make_copy(const EnsembleKalmanOutput &another);
  
  /** @brief The all ensembles. */
  std::deque<Ensemble> all_ensembles;
  
  /** @brief The times. */
  std::deque<double> times;
  
  /** @brief The llhds. */
  std::deque<double> llhds;
  
  /** @brief The lag. */
  size_t lag;
  /** @brief The output lag. */
  size_t output_lag;
  
  /** @brief The start time. */
  std::chrono::high_resolution_clock::time_point start_time;
  
};
}

#endif
