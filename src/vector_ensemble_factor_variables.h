#ifndef VECTORENSEMBLEFACTORVARIABLES_H
#define VECTORENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "ensemble_factor_variables.h"

namespace ilike
{

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

  VectorEnsembleFactorVariables(const VectorEnsembleFactorVariables &another);
  
  void operator=(const VectorEnsembleFactorVariables &another);
  EnsembleFactorVariables* duplicate() const;
  
  std::vector<arma::rowvec> get_measurement_states_for_covariance() const;
  std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const;
  std::vector<arma::colvec> get_deterministic_shifts() const;
  std::vector<arma::colvec*> get_measurements() const;
  
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double inverse_incremental_temperature,
                                             const std::vector<arma::mat> &inv_sigma_precomps,
                                             const std::vector<double> &log_det_precomps) const;

  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double inverse_incremental_temperature,
                                                       const std::vector<arma::mat> &inv_sigma_precomps,
                                                       const std::vector<double> &log_det_precomps) const;

  double evaluate_likelihoods(const Index* index) const;
  double subsample_evaluate_likelihoods(const Index* index) const;
  
  const EnsembleFactors* get_ensemble_factors() const;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;
  
  void forget_you_were_already_written_to_file();
  
  void close_ofstreams();

protected:
  
  // Stored in likelihood_estimator.
  const VectorEnsembleFactors* vector_ensemble_factors;
  
  void make_copy(const VectorEnsembleFactorVariables &another);
  
  // stored here
  std::vector<std::unique_ptr<MeasurementCovarianceEstimatorOutput>> measurement_covariance_estimator_outputs;
  
};
}

#endif
