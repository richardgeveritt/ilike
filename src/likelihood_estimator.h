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
  /**
   * @file likelihood_estimator.h
   * @brief Defines the LikelihoodEstimatorOutput class.
   *
   * Stores the output of a likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class LikelihoodEstimatorOutput
   * @brief The likelihood estimator output class.
   */


class LikelihoodEstimatorOutput;
class SMCWorker;
class Factors;

class LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the likelihoodestimator operation.
   */
  LikelihoodEstimator();
  LikelihoodEstimator(RandomNumberGenerator* rng_in,
                      size_t* seed_in,
                      Data* data_in,
                      const Parameters &algorithm_parameters_in,
                      bool smcfixed_flag_in);
  /**
   * @brief Performs the ~likelihoodestimator operation.
   */
  virtual ~LikelihoodEstimator();
  
  /**
   * @brief Performs the likelihoodestimator operation.
   *
   * @param another The LikelihoodEstimatorOutput instance to copy from.
   */
  LikelihoodEstimator(const LikelihoodEstimator &another);
  
  /**
   * @brief Assignment operator for LikelihoodEstimatorOutput.
   *
   * @param another The LikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const LikelihoodEstimator &another);
  /**
   * @brief Creates a deep copy of this LikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  virtual LikelihoodEstimator* duplicate() const=0;
  
  // Initial simulate involves simulating any of the random variables that are needed in the estimator (except for cases where we will update thesee rvs iteratively.
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  virtual LikelihoodEstimatorOutput* initialise()=0;
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  virtual LikelihoodEstimatorOutput* initialise(const Parameters &parameters)=0;
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  virtual void setup()=0;
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  virtual void setup(const Parameters &parameters)=0;
  
  // To be called if we just want the likelihood, without splitting the estimation into multiple steps.
  //double estimate();
  /**
   * @brief Performs the estimate operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double estimate(const Parameters &parameters);
  
  /**
   * @brief Performs the change data operation.
   */
  void change_data();
  /**
   * @brief Performs the change data operation.
   *
   * @param new_data The new data.
   */
  void change_data(std::shared_ptr<Data> new_data);
  /**
   * @brief Performs the change data with raw pointer operation.
   *
   * @param new_data The new data.
   */
  void change_data_with_raw_pointer(Data* new_data);
  
  /**
   * @brief Returns the data.
   *
   * @return The result.
   */
  Data* get_data() const;
  
  /**
   * @brief Returns the current data.
   *
   * @return The result.
   */
  Data* get_current_data() const;
  
  /**
   * @brief Returns the smcfixed flag.
   *
   * @return The result.
   */
  bool get_smcfixed_flag() const;
  
  Parameters algorithm_parameters;
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  virtual void specific_change_data(Data* new_data)=0;
  
  friend SMCWorker;
  
  // Not stored here. Stored in "main'.
  /** @brief The data. */
  Data* data;
  
  // not stored here
  /** @brief The current data. */
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  /** @brief The rng. */
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  /** @brief The seed. */
  size_t* seed;
  
  // Not stored here. Stored in "main'.
  /** @brief The subsampler. */
  DataSubsampler* subsampler;
  
  //ModelAndAlgorithm model_and_algorithm;
  
  // stored here
  /** @brief The factors. */
  Factors* factors;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  /** @brief The smcfixed flag. */
  bool smcfixed_flag;
  
  /**
   * @brief Copies the state of another LikelihoodEstimatorOutput into this object.
   *
   * @param another The LikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const LikelihoodEstimator &another);
  
  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};
}

#endif
