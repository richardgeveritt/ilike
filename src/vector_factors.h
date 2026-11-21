#ifndef VECTORFACTORS_H
#define VECTORFACTORS_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "factors.h"

namespace ilike
{

class LikelihoodEstimator;
class FactorVariables;
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
  bool supports_flat_generations(const Index* index) const;
  bool supports_flat_moves(const Index* index, const Index* other_index) const;
  bool supports_parallel_chunks() const;
  bool supports_parallel_copies() const;
  void redirect_random_numbers(RandomNumberGenerator* rng_in,
                               size_t* seed_in);
  size_t number_of_likelihood_estimators() const { return this->likelihood_estimators.size(); }

  // The estimator for factor i (0-based).
  LikelihoodEstimator* get_likelihood_estimator(size_t i) const;
  
  FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters) const;

  double log_maximum(const Index* index, const Parameters &parameters) const;
  void setup();
  void setup(const Parameters &conditioned_on_parameters);
  
  Data* get_current_data();
  
  // (dev/plan.md 3.21) whether any of its estimators may defer its draws,
  // found when the list is made, so that a move checks one flag
  bool can_defer() const { return this->any_defers; }
  bool advances_at_caller_target() const;

protected:

  bool any_defers = false;
  void find_deferring();
  
  void specific_change_data(Data* new_data);
  
  void make_copy(const VectorFactors &another);
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<LikelihoodEstimator>> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector< std::shared_ptr<Data> > likelihood_estimator_temp_data;
  
};
}

#endif
