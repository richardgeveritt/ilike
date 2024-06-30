#ifndef HMMFACTORS_H
#define HMMFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factors.h"
#include "parameters.h"

class ProposalKernel;
class HMMFactorVariables;
class LikelihoodEstimator;
class FactorVariables;
class Particle;

class HMMFactors : public Factors
{

public:

  HMMFactors();
  HMMFactors(ProposalKernel* transition_kernel_in,
             const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);

  virtual ~HMMFactors();

  HMMFactors(const HMMFactors &another);

  void operator=(const HMMFactors &another);
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
  
  Data* get_current_data();
  
  //std::shared_ptr<Data> get
  
protected:
  
  void specific_change_data(Data* new_data);
  
  friend HMMFactorVariables;
  // stored here
  ProposalKernel* transition_kernel;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  // stored here
  //std::vector<Data> data_time_slices;
  
  // pointers to llhd estimators
  // stored here
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  std::vector< std::shared_ptr<Data> > likelihood_estimator_temp_data;

  void make_copy(const HMMFactors &another);

};

#endif
