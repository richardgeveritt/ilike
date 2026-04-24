#ifndef HMMENSEMBLEFACTORVARIABLES_H
#define HMMENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ensemble_factor_variables.h"

namespace ilike
{
  /**
   * @file hmm_ensemble_factor_variables.h
   * @brief Defines the MeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class MeasurementCovarianceEstimatorOutput
   * @brief The measurement covariance estimator output class.
   */


class MeasurementCovarianceEstimatorOutput;
class HMMEnsembleFactors;

class HMMEnsembleFactorVariables : public EnsembleFactorVariables
{
  
public:
  
  /**
   * @brief Performs the hmmensemblefactorvariables operation.
   */
  HMMEnsembleFactorVariables();
  
  /**
   * @brief Performs the ~hmmensemblefactorvariables operation.
   */
  virtual ~HMMEnsembleFactorVariables();
  
  HMMEnsembleFactorVariables(const HMMEnsembleFactors* hmm_ensemble_factors_in,
                             const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimator_outputs_in);
  
  /**
   * @brief Performs the hmmensemblefactorvariables operation.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  HMMEnsembleFactorVariables(const HMMEnsembleFactorVariables &another);
  
  /**
   * @brief Assignment operator for MeasurementCovarianceEstimatorOutput.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const HMMEnsembleFactorVariables &another);
  /**
   * @brief Creates a deep copy of this MeasurementCovarianceEstimatorOutput object.
   *
   * @return The result.
   */
  EnsembleFactorVariables* duplicate() const;
  
  /**
   * @brief Returns the measurement states for covariance.
   *
   * @return The result.
   */
  std::vector<arma::rowvec> get_measurement_states_for_covariance() const;
  /**
   * @brief Returns the shifts.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const;
  /**
   * @brief Returns the deterministic shifts.
   *
   * @return The result.
   */
  std::vector<arma::colvec> get_deterministic_shifts() const;
  //std::vector<arma::mat> get_unconditional_measurement_covariances(const std::vector<arma::mat> &Cyys,
  //                                                                 double inverse_incremental_temperature) const;
  /**
   * @brief Returns the measurements.
   *
   * @return The result.
   */
  std::vector<arma::colvec*> get_measurements() const;
  
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double inverse_incremental_temperature,
                                             const std::vector<arma::mat> &inv_sigma_precomps,
                                             const std::vector<double> &log_det_precomps) const;
  
  /*
   double evaluate_ensemble_likelihood_ratios(const Index* index,
   double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  /*
   double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
   double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double inverse_incremental_temperature,
                                                       const std::vector<arma::mat> &inv_sigma_precomps,
                                                       const std::vector<double> &log_det_precomps) const;
  
  /**
   * @brief Evaluates the likelihoods.
   *
   * @param index The index.
   *
   * @return The result.
   */
  double evaluate_likelihoods(const Index* index) const;
  /*
   double evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  double subsample_evaluate_likelihoods(const Index* index) const;
  /*
   double subsample_evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  
  const EnsembleFactors* get_ensemble_factors() const;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /*
   void evaluate_smcfixed_part_of_measurement_covariances(const Index* index);
   void evaluate_smcfixed_part_of_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   double evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index);
   double evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   double evaluate_measurement_covariances(const Index* index);
   double evaluate_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   void subsample_evaluate_smcfixed_part_of_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   double subsample_evaluate_smcadaptive_part_given_smcfixed_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   double subsample_evaluate_measurement_covariances(const Index* index);
   double subsample_evaluate_measurement_covariances(const Index* index,
   const Parameters &conditioned_on_parameters);
   
   arma::mat direct_get_gradient_of_log(const std::string &variable);
   arma::mat direct_get_gradient_of_log(const std::string &variable,
   const Parameters &conditioned_on_parameters);
   
   arma::mat direct_subsample_get_gradient_of_log(const std::string &variable);
   arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
   const Parameters &conditioned_on_parameters);
   
   arma::mat direct_get_gradient_of_log(const Index* index,
   const std::string &variable);
   arma::mat direct_get_gradient_of_log(const Index* index,
   const std::string &variable,
   const Parameters &conditioned_on_parameters);
   
   arma::mat direct_subsample_get_gradient_of_log(const Index* index,
   const std::string &variable);
   arma::mat direct_subsample_get_gradient_of_log(const Index* index,
   const std::string &variable,
   const Parameters &conditioned_on_parameters);
   */
  
protected:
  
  // Stored in likelihood_estimator.
  /** @brief The hmm ensemble factors. */
  const HMMEnsembleFactors* hmm_ensemble_factors;
  
  /**
   * @brief Copies the state of another MeasurementCovarianceEstimatorOutput into this object.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void make_copy(const HMMEnsembleFactorVariables &another);
  
  // stored here
  /** @brief The measurement covariance estimator outputs. */
  std::vector<MeasurementCovarianceEstimatorOutput*> measurement_covariance_estimator_outputs;
  
};
}

#endif
