#ifndef VECTORENSEMBLEFACTORS_H
#define VECTORENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
#include "ensemble_factors.h"
#include "distributions.h"

class MeasurementCovarianceEstimator;
class MeasurementCovarianceEstimatorOutput;
class EnsembleFactorVariables;
class Ensemble;
//class Data;

class VectorEnsembleFactors : public EnsembleFactors
{

public:

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
                     const std::vector<arma::mat> &Cyys);
  
  std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                         const arma::mat &Ginv,
                                         const arma::mat &Ftranspose,
                                         const std::vector<arma::mat> &Vs,
                                         double inverse_incremental_temperature) const;
  
  //std::vector<arma::colvec> pack_measurements() const;
  
  double get_incremental_likelihood(Ensemble* ensemble);
  
  const EnsembleFactors* get_ensemble_factors() const;
  
  void setup();
  void setup(const Parameters &conditioned_on_parameters);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature);
  
  //void find_measurement_covariances(EnsembleKalmanOutput* simulation);
  
protected:

  void make_copy(const VectorEnsembleFactors &another);
  
  // stored here
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of arma::colvec* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector<Data*> measurement_covariance_estimator_temp_data;

};

#endif
