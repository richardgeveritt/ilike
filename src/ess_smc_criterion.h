#ifndef ESSSMCCRITERION_H
#define ESSSMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_criterion.h"

namespace ilike
{
  /**
   * @file ess_smc_criterion.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class EnsembleKalmanOutput;

class ESSSMCCriterion : public SMCCriterion
{
  
public:
  
  /**
   * @brief Performs the esssmccriterion operation.
   */
  ESSSMCCriterion();
  
  /**
   * @brief Performs the esssmccriterion operation.
   *
   * @param desired_criterion_in The desired criterion.
   */
  ESSSMCCriterion(double desired_criterion_in);
  
  /**
   * @brief Performs the ~esssmccriterion operation.
   */
  virtual ~ESSSMCCriterion();
  
  /**
   * @brief Performs the esssmccriterion operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  ESSSMCCriterion(const ESSSMCCriterion &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const ESSSMCCriterion &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  SMCCriterion* duplicate() const;
  
  /**
   * @brief Performs the operator() operation.
   *
   * @param particles The particles.
   *
   * @return The result.
   */
  double operator()(const Particles &particles) const;
  /**
   * @brief Performs the operator() operation.
   *
   * @param particles The particles.
   *
   * @return The result.
   */
  double operator()(const Ensemble &particles) const;
  
  /**
   * @brief Finds desired criterion.
   *
   * @param current_state The current state.
   */
  void find_desired_criterion(SMCOutput* current_state);
  
  /**
   * @brief Performs the subsample find desired criterion operation.
   *
   * @param current_state The current state.
   */
  void subsample_find_desired_criterion(SMCOutput* current_state);
  
  /*
   void find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters);
   void subsample_find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  /*
   void find_desired_criterion(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  bool always_positive() const;
  
protected:
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const ESSSMCCriterion &another);
  
};
}

#endif
