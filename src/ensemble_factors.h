#ifndef ENSEMBLEFACTORS_H
#define ENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file ensemble_factors.h
   * @brief Defines the Index class.
   *
   * Provides index functionality.
   *
   * @namespace ilike
   * @class Index
   * @brief The index class.
   */


class Index;
class EnsembleFactorVariables;
class Ensemble;

class EnsembleFactors
{
  
public:
  
  /**
   * @brief Performs the ensemblefactors operation.
   */
  EnsembleFactors();
  /**
   * @brief Performs the ~ensemblefactors operation.
   */
  virtual ~EnsembleFactors();
  
  /**
   * @brief Performs the ensemblefactors operation.
   *
   * @param another The Index instance to copy from.
   */
  EnsembleFactors(const EnsembleFactors &another);
  
  /**
   * @brief Assignment operator for Index.
   *
   * @param another The Index instance to copy from.
   */
  void operator=(const EnsembleFactors &another);
  /**
   * @brief Creates a deep copy of this Index object.
   *
   * @return The result.
   */
  virtual EnsembleFactors* duplicate() const=0;
  
  /**
   * @brief Sets the data.
   *
   * @param index The index.
   */
  void set_data(size_t index);
  /**
   * @brief Sets the data.
   *
   * @param index The index.
   */
  virtual void set_data(const Index* index)=0;
  
  // should be updated to return std::vector<arma::colvec>, one for each factor
  /**
   * @brief Returns the measurements.
   *
   * @return The result.
   */
  virtual std::vector<arma::colvec*> get_measurements()=0;
  
  /**
   * @brief Simulates ensemble factor variables.
   *
   * @param simulated_parameters The simulated parameters.
   *
   * @return The result.
   */
  virtual EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const=0;
  
  /*
   virtual EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const=0;
  
  /*
   virtual EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  //virtual std::vector<arma::mat> get_measurement_covariances()=0;
  //virtual std::vector<arma::mat> get_measurement_covariances(const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Performs the need cxx operation.
   *
   * @return The result.
   */
  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const std::vector<arma::mat> &Cxys,
                             const std::vector<arma::mat> &Cyys) const=0;
  
  /*
   virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
   const arma::mat &Ginv,
   const arma::mat &Ftranspose,
   const std::vector<arma::mat> &Vs,
   double inverse_incremental_temperature) const=0;
   */
  
  virtual std::vector<arma::mat> get_sqrt_adjustments(const std::vector<arma::mat> &Cxys,
                                                      const std::vector<arma::mat> &Cyys,
                                                      double inverse_incremental_temperature) const=0;
  
  virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                                 const arma::mat &Dhathalf,
                                                 const arma::mat &P,
                                                 const arma::mat &Vtranspose,
                                                 const std::vector<arma::mat> &Yhats,
                                                 double inverse_incremental_temperature) const=0;
  
  /**
   * @brief Returns the incremental likelihood.
   *
   * @param ensemble The ensemble.
   *
   * @return The result.
   */
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
  
  //void set_temperature(double temperature_in);
  //double get_temperature() const;
  //double get_inverse_incremental_temperature() const;
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  virtual void setup()=0;
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  virtual void setup(const Parameters &conditioned_on_parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                              std::vector<arma::mat> &inv_sigma_precomps,
                                              std::vector<double> &log_det_precomps) const=0;
  
protected:
  
  /** @brief The measurement names. */
  std::vector<std::string> measurement_names;
  
  // quite hacky...
  //double temperature;
  //double previous_temperature;
  
  /**
   * @brief Copies the state of another Index into this object.
   *
   * @param another The Index instance to copy from.
   */
  void make_copy(const EnsembleFactors &another);
  
};
}

#endif
