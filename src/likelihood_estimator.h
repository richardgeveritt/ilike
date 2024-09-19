#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
//#include "model_and_algorithm.h"
#include "distributions.h"
#include "data_subsampler.h"
#include "parameters.h"

namespace ilike
{
class LikelihoodEstimatorOutput;
class SMCWorker;
class Factors;

class LikelihoodEstimator
{
  
public:
  
  LikelihoodEstimator();
  LikelihoodEstimator(RandomNumberGenerator* rng_in,
                      size_t* seed_in,
                      Data* data_in,
                      const Parameters &algorithm_parameters_in,
                      bool smcfixed_flag_in);
  virtual ~LikelihoodEstimator();
  
  LikelihoodEstimator(const LikelihoodEstimator &another);
  
  void operator=(const LikelihoodEstimator &another);
  virtual LikelihoodEstimator* duplicate() const=0;
  
  // Initial simulate involves simulating any of the random variables that are needed in the estimator (except for cases where we will update thesee rvs iteratively.
  virtual LikelihoodEstimatorOutput* initialise()=0;
  virtual LikelihoodEstimatorOutput* initialise(const Parameters &parameters)=0;
  
  virtual void setup()=0;
  virtual void setup(const Parameters &parameters)=0;
  
  // To be called if we just want the likelihood, without splitting the estimation into multiple steps.
  //double estimate();
  double estimate(const Parameters &parameters);
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  void change_data_with_raw_pointer(Data* new_data);
  
  Data* get_data() const;
  
  Data* get_current_data() const;
  
  bool get_smcfixed_flag() const;
  
  Parameters algorithm_parameters;
  
protected:
  
  virtual void specific_change_data(Data* new_data)=0;
  
  friend SMCWorker;
  
  // Not stored here. Stored in "main'.
  Data* data;
  
  // not stored here
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  size_t* seed;
  
  // Not stored here. Stored in "main'.
  DataSubsampler* subsampler;
  
  //ModelAndAlgorithm model_and_algorithm;
  
  // stored here
  Factors* factors;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  void make_copy(const LikelihoodEstimator &another);
  
  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};
}

#endif
