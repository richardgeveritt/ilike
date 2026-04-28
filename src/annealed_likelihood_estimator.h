#ifndef ANNEALEDLIKELIHOODESTIMATOR_H
#define ANNEALEDLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <memory>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{
  /**
   * @file annealed_likelihood_estimator.h
   * @brief Defines the AnnealedLikelihoodEstimatorOutput class.
   *
   * Stores the output of a annealed likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class AnnealedLikelihoodEstimatorOutput
   * @brief The annealed likelihood estimator output class.
   */


class AnnealedLikelihoodEstimatorOutput;
class IndependentProposalKernel;

/*
 class Power
 {
 
 public:
 
 Power();
 ~Power();
 
 Power(const Power &another);
 
 void operator=(const Power &another);
 Power* duplicate() const;
 
 virtual double get_power() const=0;
 
 
 
 protected:
 
 void make_copy(const Power &another);
 };
 
 class FunctionPower : public Power
 {
 
 public:
 
 FunctionPower();
 ~FunctionPower();
 
 FunctionPower(const FunctionPower &another);
 
 void operator=(const FunctionPower &another);
 FunctionPower* duplicate() const;
 
 protected:
 
 void make_copy(const FunctionPower &another);
 };
 
 class ConstantPower : public Power
 {
 
 public:
 
 DoublePower();
 ~DoublePower();
 
 DoublePower(const DoublePower &another);
 
 void operator=(const DoublePower &another);
 DoublePower* duplicate() const;
 
 protected:
 
 void make_copy(const DoublePower &another);
 };
 */


class AnnealedLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the annealedlikelihoodestimator operation.
   */
  AnnealedLikelihoodEstimator();
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              PowerFunctionPtr function_power_in,
                              const std::string &power_variable_in,
                              bool smcfixed_flag_in);
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              double constant_power_in,
                              bool smcfixed_flag_in);
  
  /**
   * @brief Performs the ~annealedlikelihoodestimator operation.
   */
  virtual ~AnnealedLikelihoodEstimator();
  
  /**
   * @brief Performs the annealedlikelihoodestimator operation.
   *
   * @param another The AnnealedLikelihoodEstimatorOutput instance to copy from.
   */
  AnnealedLikelihoodEstimator(const AnnealedLikelihoodEstimator &another);
  
  /**
   * @brief Assignment operator for AnnealedLikelihoodEstimatorOutput.
   *
   * @param another The AnnealedLikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const AnnealedLikelihoodEstimator &another);
  /**
   * @brief Creates a deep copy of this AnnealedLikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend AnnealedLikelihoodEstimatorOutput;
  
  /** @brief The function power. */
  PowerFunctionPtr function_power;
  /** @brief The power variable. */
  std::string power_variable;
  /** @brief The constant power. */
  double constant_power;
  /** @brief The use constant. */
  bool use_constant;
  
  // Stored here.
  /** @brief The estimator. */
  LikelihoodEstimator* estimator;
  
  /** @brief HDF5 output file (kept open for the duration of a run). */
  std::shared_ptr<HighFive::File> h5_file;
  /** @brief Path to the HDF5 output file. */
  std::string h5_file_path;
  
  // Stored here.
  //AnnealedLikelihoodEstimatorOutput* output;
  
  /**
   * @brief Copies the state of another AnnealedLikelihoodEstimatorOutput into this object.
   *
   * @param another The AnnealedLikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const AnnealedLikelihoodEstimator &another);
  
};
}

#endif
