#ifndef ENSEMBLEFACTORS_H
#define ENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{

class Index;
class EnsembleFactorVariables;
class Ensemble;

class EnsembleFactors
{
  
public:
  
  EnsembleFactors();
  virtual ~EnsembleFactors();
  // (C270) the measurement estimators draw from this generator and seed
  // counter instead; nothing by default
  virtual void redirect_random_numbers(RandomNumberGenerator* rng_in, size_t* seed_in) {}
  // (C278) several measurement factors are assimilated one after another: the
  // number of measurement factors, and which one the next pass uses (-1: all,
  // the default for every other use)
  virtual size_t number_of_measurement_factors() const { return 1; }
  virtual void set_active_measurement_factor(int factor) {}
  // (C279) compact ensembles: whether every measurement factor supports them,
  // and factor j's predicted measurement and noise (see
  // MeasurementCovarianceEstimator)
  virtual bool supports_compact_members() const { return false; }
  // the measurement factors the update uses now (all, or the active one, C278)
  virtual std::vector<size_t> active_measurement_factor_indices() const
  {
    std::vector<size_t> all(this->number_of_measurement_factors());
    for (size_t i = 0; i < all.size(); ++i)
      all[i] = i;
    return all;
  }
  virtual arma::colvec compact_measurement_state(size_t j, const Parameters &parameters) { Rcpp::stop("EnsembleFactors::compact_measurement_state - not available."); }
  virtual arma::colvec compact_noise(size_t j) { Rcpp::stop("EnsembleFactors::compact_noise - not available."); }
  virtual void compact_simulate(size_t j, const Parameters &parameters, arma::colvec &measurement_state, arma::colvec &noise) { Rcpp::stop("EnsembleFactors::compact_simulate - not available."); }
  virtual arma::colvec compact_shift(size_t j, const arma::colvec &measurement_state, const arma::colvec &noise, double inverse_incremental_temperature) const { Rcpp::stop("EnsembleFactors::compact_shift - not available."); }
  // (C280) member i's weight from its rows: the sum over the estimators in
  // index of each one's ratio (a member's evaluate_ensemble_likelihood_ratios)
  virtual double compact_ensemble_likelihood_ratios(const Index* index,
                                                    const Ensemble* ensemble,
                                                    size_t i,
                                                    double inverse_incremental_temperature,
                                                    const std::vector<arma::mat> &inv_sigma_precomps,
                                                    const std::vector<double> &log_det_precomps) const { Rcpp::stop("EnsembleFactors::compact_ensemble_likelihood_ratios - not available."); }
  
  EnsembleFactors(const EnsembleFactors &another);
  
  void operator=(const EnsembleFactors &another);
  virtual EnsembleFactors* duplicate() const=0;
  
  void set_data(size_t index);
  virtual void set_data(const Index* index)=0;
  
  // should be updated to return std::vector<arma::colvec>, one for each factor
  virtual std::vector<arma::colvec*> get_measurements()=0;
  
  virtual EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const=0;

  virtual EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const=0;

  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const std::vector<arma::mat> &Cxys,
                             const std::vector<arma::mat> &Cyys) const=0;

  virtual std::vector<arma::mat> get_sqrt_adjustments(const std::vector<arma::mat> &Cxys,
                                                      const std::vector<arma::mat> &Cyys,
                                                      double inverse_incremental_temperature) const=0;
  
  virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                                 const arma::mat &Dhathalf,
                                                 const arma::mat &P,
                                                 const arma::mat &Vtranspose,
                                                 const std::vector<arma::mat> &Yhats,
                                                 double inverse_incremental_temperature) const=0;
  
  virtual double get_incremental_likelihood(Ensemble* ensemble) const=0;
  virtual double get_mc_inversion_incremental_likelihood(Ensemble* ensemble,
                                                      double inverse_incremental_temperature) const=0;
  virtual double get_inversion_incremental_likelihood(Ensemble* ensemble,
                                                      double inverse_incremental_temperature) const=0;
  virtual double get_unbiased_inversion_incremental_likelihood(Ensemble* ensemble,
                                                               double inverse_incremental_temperature) const=0;
  virtual void get_path1_inversion_incremental_likelihood(Ensemble* ensemble,
                                                          std::vector<double> &log_measurement_likelihood_means,
                                                          double temperature,
                                                          double multiplier) const=0;
  virtual void get_path2_inversion_incremental_likelihood(Ensemble* ensemble,
                                                          std::vector<double> &log_measurement_likelihood_means,
                                                          std::vector<double> &log_measurement_likelihood_variances) const=0;
  
  virtual void calculate_kalman_gains(Ensemble* ensemble,
                                      double inverse_incremental_temperature) const=0;

  virtual void setup()=0;
  virtual void setup(const Parameters &conditioned_on_parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                              std::vector<arma::mat> &inv_sigma_precomps,
                                              std::vector<double> &log_det_precomps) const=0;
  
protected:
  
  std::vector<std::string> measurement_names;
  
  // quite hacky...
  //double temperature;
  //double previous_temperature;
  
  void make_copy(const EnsembleFactors &another);
  
};
}

#endif
