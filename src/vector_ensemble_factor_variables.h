#ifndef VECTORENSEMBLEFACTORVARIABLES_H
#define VECTORENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ensemble_factor_variables.h"

class MeasurementCovarianceEstimatorOutput;
class EnsembleFactorVariables;
class VectorEnsembleFactors;

class VectorEnsembleFactorVariables : public EnsembleFactorVariables
{

public:

  VectorEnsembleFactorVariables();
  
  VectorEnsembleFactorVariables(const VectorEnsembleFactors* ensemble_factors_in,
                                const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimator_outputs_in);

  virtual ~VectorEnsembleFactorVariables();
  
  //VectorEnsembleFactorVariables(std::vector<MeasurementCovarianceEstimatorOutput*> measurement_covariance_estimator_outputs_in);

  VectorEnsembleFactorVariables(const VectorEnsembleFactorVariables &another);

  void operator=(const VectorEnsembleFactorVariables &another);
  EnsembleFactorVariables* duplicate() const;
  
  std::vector<arma::rowvec> get_measurement_states_for_covariance() const;
  std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const;
  std::vector<arma::colvec> get_deterministic_shifts() const;
  //std::vector<arma::mat> get_unconditional_measurement_covariances(const std::vector<arma::mat> &Cyys,
  //                                                                 double inverse_incremental_temperature) const;
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
  
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double inverse_incremental_temperature,
                                                       const std::vector<arma::mat> &inv_sigma_precomps,
                                                       const std::vector<double> &log_det_precomps) const;
  
  /*
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double inverse_incremental_temperature,
                                                       const Parameters &conditioned_on_parameters);
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
  
  void forget_you_were_already_written_to_file();
  
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
  const VectorEnsembleFactors* vector_ensemble_factors;

  void make_copy(const VectorEnsembleFactorVariables &another);
  
  // stored here
  std::vector<MeasurementCovarianceEstimatorOutput*> measurement_covariance_estimator_outputs;

};

#endif
