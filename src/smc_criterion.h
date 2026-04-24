#ifndef SMCCRITERION_H
#define SMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"
#include "ensemble.h"

namespace ilike
{
  /**
   * @file smc_criterion.h
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

class SMCCriterion
{
  
public:
  
  /**
   * @brief Performs the smccriterion operation.
   */
  SMCCriterion();
  /**
   * @brief Performs the ~smccriterion operation.
   */
  virtual ~SMCCriterion();
  
  /**
   * @brief Performs the smccriterion operation.
   *
   * @param desired_criterion_in The desired criterion.
   */
  SMCCriterion(double desired_criterion_in);
  
  /**
   * @brief Performs the smccriterion operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  SMCCriterion(const SMCCriterion &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const SMCCriterion &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  virtual SMCCriterion* duplicate() const=0;
  
  /**
   * @brief Performs the operator() operation.
   *
   * @param particles The particles.
   *
   * @return The result.
   */
  virtual double operator()(const Particles &particles) const=0;
  
  /**
   * @brief Finds desired criterion.
   *
   * @param current_state The current state.
   */
  virtual void find_desired_criterion(SMCOutput* current_state)=0;
  
  /**
   * @brief Performs the subsample find desired criterion operation.
   *
   * @param current_state The current state.
   */
  virtual void subsample_find_desired_criterion(SMCOutput* current_state)=0;
  
  /*
   virtual void find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual double operator()(const Ensemble &particles) const=0;
  
  /**
   * @brief Finds desired criterion.
   *
   * @param current_state The current state.
   */
  virtual void find_desired_criterion(EnsembleKalmanOutput* current_state)=0;
  
  /*
   virtual void find_desired_criterion(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual bool always_positive() const=0;
  
protected:
  
  /** @brief The desired criterion. */
  double desired_criterion;
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const SMCCriterion &another);
  
};
}

#endif
