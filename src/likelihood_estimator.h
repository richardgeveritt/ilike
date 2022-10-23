#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
//#include "model_and_algorithm.h"
#include "distributions.h"
#include "data_subsampler.h"
#include "data.h"

class LikelihoodEstimatorOutput;
class SMCWorker;
class Factors;

class LikelihoodEstimator
{

public:

  LikelihoodEstimator();
  LikelihoodEstimator(RandomNumberGenerator* rng_in,
                      size_t* seed_in,
                      Data* data_in);
  virtual ~LikelihoodEstimator();

  LikelihoodEstimator(const LikelihoodEstimator &another);

  void operator=(const LikelihoodEstimator &another);
  virtual LikelihoodEstimator* duplicate() const=0;

  // Initial simulate involves simulating any of the random variables that are needed in the estimator (except for cases where we will update thesee rvs iteratively.
  virtual LikelihoodEstimatorOutput* initialise()=0;
  virtual LikelihoodEstimatorOutput* initialise(const Parameters &parameters)=0;
  
  // To be called if we just want the likelihood, without splitting the estimation into multiple steps.
  //double estimate();
  double estimate(const Parameters &parameters);
  
  void change_data();
  void change_data(Data* new_data);
  
  Data* get_data() const;

protected:

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

  void make_copy(const LikelihoodEstimator &another);

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;

  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
