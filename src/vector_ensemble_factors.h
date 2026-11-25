#ifndef VECTORENSEMBLEFACTORS_H
#define VECTORENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
#include "ensemble_factors.h"
#include "distributions.h"

namespace ilike
{

class MeasurementCovarianceEstimator;
class MeasurementCovarianceEstimatorOutput;
class EnsembleFactorVariables;
class Ensemble;

class VectorEnsembleFactors : public EnsembleFactors
{
  
public:
  // (C280) compact ensembles
  void compact_simulate(size_t j, const Parameters &parameters, arma::colvec &measurement_state, arma::colvec &noise);
  arma::colvec compact_shift(size_t j, const arma::colvec &measurement_state, const arma::colvec &noise, double inverse_incremental_temperature) const;
  double compact_ensemble_likelihood_ratios(const Index* index,
                                            const Ensemble* ensemble,
                                            size_t i,
                                            double inverse_incremental_temperature,
                                            const std::vector<arma::mat> &inv_sigma_precomps,
                                            const std::vector<double> &log_det_precomps) const;
  size_t number_of_measurement_factors() const { return this->measurement_covariance_estimators.size(); }
  bool supports_compact_members() const;
  arma::colvec compact_measurement_state(size_t j, const Parameters &parameters);
  arma::colvec compact_noise(size_t j);

  
  VectorEnsembleFactors();
  VectorEnsembleFactors(const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in);
  
  virtual ~VectorEnsembleFactors();
  
  VectorEnsembleFactors(const VectorEnsembleFactors &another);
  
  void operator=(const VectorEnsembleFactors &another);
  EnsembleFactors* duplicate() const;
  
  void set_data(const Index* index);
  
  // should be updated to return std::vector<arma::colvec>, one for each factor
  std::vector<arma::colvec*> get_measurements();
  
  EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const;
  EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const;

  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const std::vector<arma::mat> &Cxys,
                     const std::vector<arma::mat> &Cyys) const;
  
  std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                         const arma::mat &Dhathalf,
                                         const arma::mat &P,
                                         const arma::mat &Vtranspose,
                                         const std::vector<arma::mat> &Yhats,
                                         double inverse_incremental_temperature) const;
  
  std::vector<arma::mat> get_sqrt_adjustments(const std::vector<arma::mat> &Cxys,
                                              const std::vector<arma::mat> &Cyys,
                                              double inverse_incremental_temperature) const;

  double get_incremental_likelihood(Ensemble* ensemble) const;
  double get_mc_inversion_incremental_likelihood(Ensemble* ensemble,
                                              double inverse_incremental_temperature) const;
  double get_inversion_incremental_likelihood(Ensemble* ensemble,
                                              double inverse_incremental_temperature) const;
  double get_unbiased_inversion_incremental_likelihood(Ensemble* ensemble,
                                                       double inverse_incremental_temperature) const;
  void get_path1_inversion_incremental_likelihood(Ensemble* ensemble,
                                                  std::vector<double> &log_measurement_likelihood_means,
                                                  double temperature,
                                                  double multiplier) const;
  void get_path2_inversion_incremental_likelihood(Ensemble* ensemble,
                                                  std::vector<double> &log_measurement_likelihood_means,
                                                  std::vector<double> &log_measurement_likelihood_variances) const;
  
  void calculate_kalman_gains(Ensemble* ensemble,
                              double inverse_incremental_temperature) const;
  
  const EnsembleFactors* get_ensemble_factors() const;
  
  void setup();
  void setup(const Parameters &conditioned_on_parameters);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                      std::vector<arma::mat> &inv_sigma_precomps,
                                      std::vector<double> &log_det_precomps) const;

protected:
  
  void make_copy(const VectorEnsembleFactors &another);
  
  // stored here
  std::vector<std::unique_ptr<MeasurementCovarianceEstimator>> measurement_covariance_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of arma::colvec* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector< std::shared_ptr<Data> > measurement_covariance_estimator_temp_data;
  
};
}

#endif
