#ifndef ENSEMBLEKALMANOUTPUT_H
#define ENSEMBLEKALMANOUTPUT_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"
#include "ensemble.h"

class EnsembleKalman;
class EnsembleKalmanFilter;
class EnsembleSequencer;
class IterativeEnsembleKalmanInversion;

class EnsembleKalmanOutput : public LikelihoodEstimatorOutput
{

public:

  // unsure
  // copied from KF
  
  EnsembleKalmanOutput();
  EnsembleKalmanOutput(EnsembleKalman* estimator_in,
                       size_t lag_in);
  virtual ~EnsembleKalmanOutput();

  EnsembleKalmanOutput(const EnsembleKalmanOutput &another);
  void operator=(const EnsembleKalmanOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  Ensemble* add_ensemble();
  //void add_proposed_ensemble(const Ensemble &latest_proposals);
  //void initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  void initialise_next_step();
  
  Ensemble back() const;
  Ensemble& back();
  
  void simulate();
  
  void simulate(const Parameters &parameters);
  
  //double evaluate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters);
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  size_t number_of_ensemble_kalman_iterations() const;
  
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

protected:
  
  friend EnsembleKalmanFilter;
  friend IterativeEnsembleKalmanInversion;
  friend EnsembleSequencer;

  // Stored in Factors.
  EnsembleKalman* estimator;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;

  void make_copy(const EnsembleKalmanOutput &another);
  
  std::deque<Ensemble> all_ensembles;
  
  size_t lag;

};

#endif