#ifndef HMMENSEMBLEFACTORS_H
#define HMMENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "ensemble_factors.h"

namespace ilike
{

class MeasurementCovarianceEstimator;
class EnsembleFactorVariables;
class ProposalKernel;
class Ensemble;

class HMMEnsembleFactors : public EnsembleFactors
{
  friend class HMMEnsembleFactorVariables; // (C265: repointing to a copy's estimators)
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

  void redirect_random_numbers(RandomNumberGenerator* rng_in, size_t* seed_in); // (C270)
  // (C278) see EnsembleFactors; the loops of the filter's update (gains,
  // likelihood, adjustments, measurements, and the members' measurement
  // states and shifts) take only the active factor's estimator
  size_t number_of_measurement_factors() const { return this->measurement_covariance_estimators.size(); }
  void set_active_measurement_factor(int factor) { this->active_measurement_factor = factor; }
  const std::vector<size_t>& active_measurement_factors() const;
  bool supports_compact_members() const; // (C279)
  std::vector<size_t> active_measurement_factor_indices() const { return this->active_measurement_factors(); }
  arma::colvec compact_measurement_state(size_t j, const Parameters &parameters);
  arma::colvec compact_noise(size_t j);
private:
  int active_measurement_factor = -1;
  // (the list active_measurement_factors returns, kept so that the loops per
  // member allocate nothing)
  mutable std::vector<size_t> active_factors;
  
public:
  
  HMMEnsembleFactors();
  
  HMMEnsembleFactors(ProposalKernel* transition_kernel_in,
                     const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in);
  
  virtual ~HMMEnsembleFactors();
  
  HMMEnsembleFactors(const HMMEnsembleFactors &another);
  
  void operator=(const HMMEnsembleFactors &another);
  EnsembleFactors* duplicate() const;

  ProposalKernel* get_transition_kernel() const;
  
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
  
  void setup();
  void setup(const Parameters &conditioned_on_parameters);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                      std::vector<arma::mat> &inv_sigma_precomps,
                                      std::vector<double> &log_det_precomps) const;
  
protected:
  
  void make_copy(const HMMEnsembleFactors &another);
  
  // stored here
  std::unique_ptr<ProposalKernel> transition_kernel;
  
  // stored here
  std::vector<std::unique_ptr<MeasurementCovarianceEstimator>> measurement_covariance_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector< std::shared_ptr<Data> > measurement_covariance_estimator_temp_data;
  
};
}

#endif
