#ifndef HMMENSEMBLEFACTORVARIABLES_H
#define HMMENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ensemble_factor_variables.h"

class MeasurementCovarianceEstimatorOutput;
class HMMEnsembleFactors;

class HMMEnsembleFactorVariables : public EnsembleFactorVariables
{

public:

  HMMEnsembleFactorVariables();

  virtual ~HMMEnsembleFactorVariables();
  
  HMMEnsembleFactorVariables(HMMEnsembleFactors* hmm_ensemble_factors_in,
                             const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimator_outputs_in);

  HMMEnsembleFactorVariables(const HMMEnsembleFactorVariables &another);

  void operator=(const HMMEnsembleFactorVariables &another);
  EnsembleFactorVariables* duplicate() const;
  
  std::vector<arma::rowvec> get_measurement_states_for_covariance() const;
  std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const;
  std::vector<arma::colvec> get_deterministic_shifts() const;
  std::vector<arma::mat> get_kalman_gains(const std::vector<arma::mat> &Cxys,
                                          const std::vector<arma::mat> &Cyys,
                                          double inverse_incremental_temperature) const;
  std::vector<arma::colvec*> get_measurements() const;
  
  std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                         const arma::mat &Ginv,
                                         const arma::mat &Ftranspose,
                                         const std::vector<arma::mat> &Vs,
                                         double inverse_incremental_temperature) const;
  
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double inverse_incremental_temperature);
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double inverse_incremental_temperature,
                                             const Parameters &conditioned_on_parameters);
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double inverse_incremental_temperature,
                                                       const Parameters &conditioned_on_parameters);
  
  double evaluate_likelihoods(const Index* index);
  double evaluate_likelihoods(const Index* index,
                              const Parameters &conditioned_on_parameters);
  double subsample_evaluate_likelihoods(const Index* index);
  double subsample_evaluate_likelihoods(const Index* index,
                                        const Parameters &conditioned_on_parameters);
  
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

  void make_copy(const HMMEnsembleFactorVariables &another);
  
  // stored here
  std::vector<MeasurementCovarianceEstimatorOutput*> measurement_covariance_estimator_outputs;

};

#endif
