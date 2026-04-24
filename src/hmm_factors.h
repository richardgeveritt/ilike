#ifndef HMMFACTORS_H
#define HMMFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factors.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file hmm_factors.h
   * @brief Defines the ProposalKernel class.
   *
   * A generic proposal kernel. Proposes new parameter values during MCMC or SMC moves using a generic distribution centred on the current state.
   *
   * @namespace ilike
   * @class ProposalKernel
   * @brief The proposal kernel class.
   */


class ProposalKernel;
class HMMFactorVariables;
class LikelihoodEstimator;
class FactorVariables;
class Particle;

class HMMFactors : public Factors
{
  
public:
  
  /**
   * @brief Performs the hmmfactors operation.
   */
  HMMFactors();
  HMMFactors(ProposalKernel* transition_kernel_in,
             const std::vector<LikelihoodEstimator*> &likelihood_estimators_in);
  
  /**
   * @brief Performs the ~hmmfactors operation.
   */
  virtual ~HMMFactors();
  
  /**
   * @brief Performs the hmmfactors operation.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  HMMFactors(const HMMFactors &another);
  
  /**
   * @brief Assignment operator for ProposalKernel.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void operator=(const HMMFactors &another);
  /**
   * @brief Creates a deep copy of this ProposalKernel object.
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
  
  //std::shared_ptr<Data> get
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend HMMFactorVariables;
  // stored here
  /** @brief The transition kernel. */
  ProposalKernel* transition_kernel;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  /** @brief The smcfixed flag. */
  bool smcfixed_flag;
  
  // stored here
  //std::vector<Data> data_time_slices;
  
  // pointers to llhd estimators
  // stored here
  /** @brief The likelihood estimators. */
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // stored here
  // data temporarily used in a likelihood estimator
  /** @brief The likelihood estimator temp data. */
  std::vector< std::shared_ptr<Data> > likelihood_estimator_temp_data;
  
  /**
   * @brief Copies the state of another ProposalKernel into this object.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void make_copy(const HMMFactors &another);
  
};
}

#endif
