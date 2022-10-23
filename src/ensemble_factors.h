#ifndef ENSEMBLEFACTORS_H
#define ENSEMBLEFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

class Index;
class EnsembleFactorVariables;

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
  virtual arma::colvec get_measurements()=0;
  
  virtual EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters)=0;
  virtual EnsembleFactorVariables* simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
                                                                      const Parameters &conditioned_on_parameters)=0;
  
  //virtual std::vector<arma::mat> get_measurement_covariances()=0;
  //virtual std::vector<arma::mat> get_measurement_covariances(const Parameters &conditioned_on_parameters)=0;

  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const std::vector<arma::mat> &Cxys,
                             const std::vector<arma::mat> &Cyys)=0;
  
  void set_temperature(double temperature_in);
  double get_temperature() const;
  
protected:

  std::vector<std::string> measurement_names;
  
  // quite hacky...
  double temperature;
  
  void make_copy(const EnsembleFactors &another);

};

#endif
