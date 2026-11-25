#ifndef KALMANFILTEROUTPUT_H
#define KALMANFILTEROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

class KalmanFilter;

class KalmanFilterOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  KalmanFilterOutput();
  KalmanFilterOutput(KalmanFilter* estimator_in,
                     size_t lag_in,
                     const std::string &results_name_in);
  virtual ~KalmanFilterOutput();
  
  KalmanFilterOutput(const KalmanFilterOutput &another);
  void operator=(const KalmanFilterOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  
  void simulate();
  
  void simulate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters);
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void set_current_predicted_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  void add_predicted_statistics();
  
  void set_current_posterior_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  void add_posterior_statistics();
  
  void add_schedule_parameters(const Parameters &current_schedule_parameters);
  
  void add_log_normalising_constant_ratio(double log_normalising_constant_ratio);
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  void set_time();
  void set_llhd(double llhd_in);
  
  arma::colvec predicted_mean_back() const;
  arma::colvec posterior_mean_back() const;
  arma::mat predicted_covariance_back() const;
  arma::mat posterior_covariance_back() const;
  
  void set_current_predicted_to_be_current_posterior();
  
  size_t predicted_size() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  void forget_you_were_already_written_to_file();
  
  void terminate();
  
  void increment_kf_iteration();
  // (dev/plan.md 3.21 (c2''')) following the caller's data tempering: the
  // observations committed, their log-likelihood, and whether this output has
  // been run through the caller's committed observations yet
  size_t followed_steps = 0;
  double committed_log_likelihood = 0.0;
  bool caught_up = false;
  void advance_at_caller_target(const Parameters &parameters);
  
  // (C264) a compact record when nothing is stored (iterations_to_store =
  // 0, the default inside another algorithm): the log-likelihood, and the
  // file write_to_file would make (its attributes only); an MCMC chain then
  // stores each iteration as numbers instead of a copy of the particle
  size_t compact_size() const { return (this->output_lag == 0) ? 1 : 0; }
  void compact_values(double *values) const { values[0] = this->log_likelihood; }
  void write_compact_values(const std::string &directory_name,
                            const std::string &index,
                            const double *values);
  void close_ofstreams();
  
  void print(std::ostream &os) const;
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
  std::string results_name;
  
  std::chrono::high_resolution_clock::time_point start_time;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  KalmanFilter* estimator;
  HighFive::File& output_file(const std::string &dir_name);
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index="");
  
  void make_copy(const KalmanFilterOutput &another);
  
  std::deque<arma::colvec> predicted_means;
  std::deque<arma::mat> predicted_covariances;
  std::deque<arma::colvec> posterior_means;
  std::deque<arma::mat> posterior_covariances;
  
  std::deque<Parameters> schedule_parameters;
  
  std::deque<double> log_normalising_constant_ratios;
  
  arma::colvec current_predicted_mean;
  arma::mat current_predicted_covariance;
  arma::colvec current_posterior_mean;
  arma::mat current_posterior_covariance;
  
  size_t lag;
  size_t output_lag;
  
  size_t kf_iteration;
  
  int iteration_written_to_file;
  
};
}

#endif
