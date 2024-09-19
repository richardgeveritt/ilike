#ifndef HMMENSEMBLEFACTORS_H
#define HMMENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ensemble_factors.h"

namespace ilike
{
class MeasurementCovarianceEstimator;
class EnsembleFactorVariables;
class ProposalKernel;
class Ensemble;
//class Data;

class HMMEnsembleFactors : public EnsembleFactors
{
  
public:
  
  HMMEnsembleFactors();
  
  HMMEnsembleFactors(ProposalKernel* transition_kernel_in,
                     const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in);
  
  virtual ~HMMEnsembleFactors();
  
  HMMEnsembleFactors(const HMMEnsembleFactors &another);
  
  void operator=(const HMMEnsembleFactors &another);
  EnsembleFactors* duplicate() const;
  
  void set_data(const Index* index);
  
  // should be updated to return std::vector<arma::colvec>, one for each factor
  std::vector<arma::colvec*> get_measurements();
  
  EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const;
  /*
   EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
   const Parameters &conditioned_on_parameters);
   */
  EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const;
  /*
   EnsembleFactorVariables* subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
   const Parameters &conditioned_on_parameters);
   */
  
  //std::vector<arma::mat> get_measurement_covariances();
  //std::vector<arma::mat> get_measurement_covariances(const Parameters &conditioned_on_parameters);
  
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
  ProposalKernel* transition_kernel;
  
  // stored here
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector< std::shared_ptr<Data> > measurement_covariance_estimator_temp_data;
  
};
}

#endif
