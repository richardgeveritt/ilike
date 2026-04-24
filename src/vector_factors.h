#ifndef VECTORFACTORS_H
#define VECTORFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factors.h"

namespace ilike
{
  /**
   * @file vector_factors.h
   * @brief Defines the LikelihoodEstimator class.
   *
   * Estimates the likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class LikelihoodEstimator
   * @brief The likelihood estimator class.
   */


class LikelihoodEstimator;
class FactorVariables;
//class Data;
class Particle;

class VectorFactors : public Factors
{
  
public:
  
  /**
   * @brief Performs the vectorfactors operation.
   */
  VectorFactors();
  /**
   * @brief Performs the vectorfactors operation.
   *
   * @param likelihood_estimators_in The likelihood estimators.
   */
  VectorFactors(const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);
  
  /**
   * @brief Performs the ~vectorfactors operation.
   */
  virtual ~VectorFactors();
  
  /**
   * @brief Performs the vectorfactors operation.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  VectorFactors(const VectorFactors &another);
  
  /**
   * @brief Assignment operator for LikelihoodEstimator.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  void operator=(const VectorFactors &another);
  /**
   * @brief Creates a deep copy of this LikelihoodEstimator object.
   *
   * @return The result.
   */
  Factors* duplicate() const;
  
  /**
   * @brief Sets the data.
   *
   * @param index The index.
   */
  void set_data(const Index* index);
  
  /**
   * @brief Simulates factor variables.
   *
   * @param simulated_parameters The simulated parameters.
   *
   * @return The result.
   */
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
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void setup(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Returns the current data.
   *
   * @return The result.
   */
  Data* get_current_data();
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  /**
   * @brief Copies the state of another LikelihoodEstimator into this object.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  void make_copy(const VectorFactors &another);
  
  // stored here
  /** @brief The likelihood estimators. */
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  /** @brief The likelihood estimator temp data. */
  std::vector< std::shared_ptr<Data> > likelihood_estimator_temp_data;
  
};
}

#endif
