#ifndef VECTORFACTORS_H
#define VECTORFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factors.h"

class LikelihoodEstimator;
class FactorVariables;
//class Data;
class Particle;

class VectorFactors : public Factors
{

public:

  VectorFactors();
  VectorFactors(const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);

  virtual ~VectorFactors();

  VectorFactors(const VectorFactors &another);

  void operator=(const VectorFactors &another);
  Factors* duplicate() const;
  
  void set_data(const Index* index);
  
  FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters) const;
  /*
  FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters,
                                             const Parameters &conditioned_on_parameters);
  */
   
  FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters) const;
   
  /*
  FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters,
                                                       const Parameters &conditioned_on_parameters);
  */
  
  void setup();
  void setup(const Parameters &conditioned_on_parameters);
  
protected:

  void make_copy(const VectorFactors &another);
  
  // stored here
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector< std::shared_ptr<Data> > likelihood_estimator_temp_data;

};

#endif
