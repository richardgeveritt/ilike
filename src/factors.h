#ifndef FACTORS_H
#define FACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{

class Index;
class FactorVariables;

class Factors
{
  
public:
  
  Factors();
  virtual ~Factors();
  
  Factors(const Factors &another);
  
  void operator=(const Factors &another);
  virtual Factors* duplicate() const=0;
  
  void set_data(size_t index);
  virtual void set_data(const Index* index)=0;
  
  void change_data(Data* new_data);
  
  virtual Data* get_current_data()=0;
  
  virtual FactorVariables* simulate_factor_variables(const Parameters &simulated_parameters) const=0;

  // True if every estimator of the index supports flat generations
  // (LikelihoodEstimator::supports_flat_generations). Default: false.
  virtual bool supports_flat_generations(const Index* index) const { return false; }

  // True if every factor is in index or in other_index (other_index may be
  // NULL), the estimators of index support restore_from_summary and those of
  // other_index support reinitialise (dev/plan.md 2.6a, G3). Default: false.
  virtual bool supports_flat_moves(const Index* index, const Index* other_index) const { return false; }

  // True if every estimator supports_parallel_chunks (G4). Default: false.
  virtual bool supports_parallel_chunks() const { return false; }

  // G5: every estimator supports_parallel_copies; redirect every estimator's
  // random numbers. Defaults: no; nothing.
  virtual bool supports_parallel_copies() const { return false; }
  virtual void redirect_random_numbers(RandomNumberGenerator* rng_in,
                                       size_t* seed_in) {}

  virtual FactorVariables* subsample_simulate_factor_variables(const Parameters &simulated_parameters) const=0;

  virtual void setup()=0;
  virtual void setup(const Parameters &conditioned_on_parameters)=0;
  
protected:
  
  virtual void specific_change_data(Data* new_data)=0;
  
  void make_copy(const Factors &another);
  
};
}

#endif
