#ifndef ENSEMBLEFACTORS_H
#define ENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

class Index;
class EnsembleFactorVariables;
class Ensemble;

class EnsembleFactors
{

public:

  EnsembleFactors();
  virtual ~EnsembleFactors();
  
  EnsembleFactors(const EnsembleFactors &another);

  void operator=(const EnsembleFactors &another);
  virtual EnsembleFactors* duplicate() const=0;
  
  void set_data(size_t index);
  virtual void set_data(const Index* index)=0;
  
  // should be updated to return std::vector<arma::colvec>, one for each factor
  virtual std::vector<arma::colvec*> get_measurements()=0;
  
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

  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const std::vector<arma::mat> &Cxys,
                             const std::vector<arma::mat> &Cyys)=0;
  
  /*
  virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                                 const arma::mat &Ginv,
                                                 const arma::mat &Ftranspose,
                                                 const std::vector<arma::mat> &Vs,
                                                 double inverse_incremental_temperature) const=0;
  */
  
  virtual std::vector<arma::mat> get_sqrt_adjustments(const arma::mat &Sigma,
                                                      const std::vector<arma::mat> &HSigmaHts,
                                                      double inverse_incremental_temperature) const=0;
  
  virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                                 const arma::mat &Dhathalf,
                                                 const arma::mat &P,
                                                 const arma::mat &Vtranspose,
                                                 const std::vector<arma::mat> &Yhats,
                                                 double inverse_incremental_temperature) const=0;
  
  virtual double get_incremental_likelihood(Ensemble* ensemble)=0;
  virtual double get_inversion_incremental_likelihood(Ensemble* ensemble,
                                                      double inverse_incremental_temperature)=0;
  
  //void set_temperature(double temperature_in);
  //double get_temperature() const;
  //double get_inverse_incremental_temperature() const;
  
  virtual void setup()=0;
  virtual void setup(const Parameters &conditioned_on_parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature)=0;
  
protected:

  std::vector<std::string> measurement_names;
  
  // quite hacky...
  //double temperature;
  //double previous_temperature;
  
  void make_copy(const EnsembleFactors &another);

};

#endif
