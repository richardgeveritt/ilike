#ifndef HMMFACTORS_H
#define HMMFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factors.h"
#include "data.h"

class ProposalKernel;
class HMMFactorVariables;
class LikelihoodEstimator;
class FactorVariables;
class Particle;

class HMMFactors : public Factors
{

public:

  HMMFactors();
  HMMFactors(ProposalKernel* transition_kernel_in);

  virtual ~HMMFactors();

  HMMFactors(const HMMFactors &another);

  void operator=(const HMMFactors &another);
  Factors* duplicate() const;
  
  void set_data(const Index* index);
  
  FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters);
  FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters,
                                             const Parameters &conditioned_on_parameters);
  FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters);
  FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters,
                                                       const Parameters &conditioned_on_parameters);
  
protected:
  
  friend HMMFactorVariables;
  // stored here
  ProposalKernel* transition_kernel;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  // stored here
  std::vector<Data> data_time_slices;
  
  // pointers to llhd estimators
  // stored here
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  std::vector<Data*> likelihood_estimator_temp_data;

  void make_copy(const HMMFactors &another);

};

#endif