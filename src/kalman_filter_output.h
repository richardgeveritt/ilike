#ifndef KALMANFILTEROUTPUT_H
#define KALMANFILTEROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"

namespace ilike
{
  /**
   * @file kalman_filter_output.h
   * @brief Defines the KalmanFilter class.
   *
   * Implements kalman filter. Performs Gaussian state-space inference using Kalman recursions.
   *
   * @namespace ilike
   * @class KalmanFilter
   * @brief The kalman filter class.
   */


class KalmanFilter;

class KalmanFilterOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the kalmanfilteroutput operation.
   */
  KalmanFilterOutput();
  KalmanFilterOutput(KalmanFilter* estimator_in,
                     size_t lag_in,
                     const std::string &results_name_in);
  /**
   * @brief Performs the ~kalmanfilteroutput operation.
   */
  virtual ~KalmanFilterOutput();
  
  /**
   * @brief Performs the kalmanfilteroutput operation.
   *
   * @param another The KalmanFilter instance to copy from.
   */
  KalmanFilterOutput(const KalmanFilterOutput &another);
  /**
   * @brief Assignment operator for KalmanFilter.
   *
   * @param another The KalmanFilter instance to copy from.
   */
  void operator=(const KalmanFilterOutput &another);
  /**
   * @brief Creates a deep copy of this KalmanFilter object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* duplicate() const;
  
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
  /**
   * @brief Evaluates the smcfixed part.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
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
  
  void set_current_predicted_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  /**
   * @brief Performs the add predicted statistics operation.
   */
  void add_predicted_statistics();
  
  void set_current_posterior_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  /**
   * @brief Performs the add posterior statistics operation.
   */
  void add_posterior_statistics();
  
  /**
   * @brief Performs the add schedule parameters operation.
   *
   * @param current_schedule_parameters The current schedule parameters.
   */
  void add_schedule_parameters(const Parameters &current_schedule_parameters);
  
  /**
   * @brief Performs the add log normalising constant ratio operation.
   *
   * @param log_normalising_constant_ratio The log normalising constant ratio.
   */
  void add_log_normalising_constant_ratio(double log_normalising_constant_ratio);
  
  /**
   * @brief Returns the likelihood estimator.
   *
   * @return The result.
   */
  LikelihoodEstimator* get_likelihood_estimator() const;
  
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
  
  /**
   * @brief Performs the predicted mean back operation.
   *
   * @return The result.
   */
  arma::colvec predicted_mean_back() const;
  /**
   * @brief Performs the posterior mean back operation.
   *
   * @return The result.
   */
  arma::colvec posterior_mean_back() const;
  /**
   * @brief Performs the predicted covariance back operation.
   *
   * @return The result.
   */
  arma::mat predicted_covariance_back() const;
  /**
   * @brief Performs the posterior covariance back operation.
   *
   * @return The result.
   */
  arma::mat posterior_covariance_back() const;
  
  /**
   * @brief Sets the current predicted to be current posterior.
   */
  void set_current_predicted_to_be_current_posterior();
  
  /**
   * @brief Performs the predicted size operation.
   *
   * @return The result.
   */
  size_t predicted_size() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   */
  void terminate();
  
  /**
   * @brief Performs the increment kf iteration operation.
   */
  void increment_kf_iteration();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /**
   * @brief Prints the object's state to an output stream.
   *
   * @param os The os.
   */
  void print(std::ostream &os) const;
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
  std::string results_name;
  
  std::chrono::high_resolution_clock::time_point start_time;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  /** @brief The estimator. */
  KalmanFilter* estimator;
  
  /** @brief The log likelihood smcfixed part. */
  double log_likelihood_smcfixed_part;
  /** @brief The subsample log likelihood smcfixed part. */
  double subsample_log_likelihood_smcfixed_part;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index="");
  
  /**
   * @brief Copies the state of another KalmanFilter into this object.
   *
   * @param another The KalmanFilter instance to copy from.
   */
  void make_copy(const KalmanFilterOutput &another);
  
  /** @brief The predicted means. */
  std::deque<arma::colvec> predicted_means;
  /** @brief The predicted covariances. */
  std::deque<arma::mat> predicted_covariances;
  /** @brief The posterior means. */
  std::deque<arma::colvec> posterior_means;
  /** @brief The posterior covariances. */
  std::deque<arma::mat> posterior_covariances;
  
  /** @brief The schedule parameters. */
  std::deque<Parameters> schedule_parameters;
  
  /** @brief The log normalising constant ratios. */
  std::deque<double> log_normalising_constant_ratios;
  
  /** @brief The current predicted mean. */
  arma::colvec current_predicted_mean;
  /** @brief The current predicted covariance. */
  arma::mat current_predicted_covariance;
  /** @brief The current posterior mean. */
  arma::colvec current_posterior_mean;
  /** @brief The current posterior covariance. */
  arma::mat current_posterior_covariance;
  
  /** @brief The lag. */
  size_t lag;
  /** @brief The output lag. */
  size_t output_lag;
  
  /** @brief The kf iteration. */
  size_t kf_iteration;
  
  /** @brief The iteration written to file. */
  int iteration_written_to_file;
  
};
}

#endif
