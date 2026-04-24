#ifndef POSITIVESMCCRITERION_H
#define POSITIVESMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_criterion.h"

namespace ilike
{
  /**
   * @file positive_smc_criterion.h
   * @brief Defines the PositiveSMCCriterion class.
   *
   * Implements an SMC resampling/adaptation criterion based on positive. Used by the SMC algorithm to decide when to resample particles.
   *
   * @namespace ilike
   * @class PositiveSMCCriterion
   * @brief A positive smc criterion derived from SMCCriterion.
   */


class PositiveSMCCriterion : public SMCCriterion
{
  
public:
  
  /**
   * @brief Default constructor for PositiveSMCCriterion.
   */
  PositiveSMCCriterion();
  
  /**
   * @brief Destructor for PositiveSMCCriterion.
   */
  virtual ~PositiveSMCCriterion();
  
  /**
   * @brief Copy constructor for PositiveSMCCriterion.
   *
   * @param another The PositiveSMCCriterion instance to copy from.
   */
  PositiveSMCCriterion(const PositiveSMCCriterion &another);
  
  /**
   * @brief Assignment operator for PositiveSMCCriterion.
   *
   * @param another The PositiveSMCCriterion instance to copy from.
   */
  void operator=(const PositiveSMCCriterion &another);
  /**
   * @brief Creates a deep copy of this PositiveSMCCriterion object.
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
  
  double operator()(const Ensemble &particles) const;
  
  /**
   * @brief Finds desired criterion.
   *
   * @param current_state The current state.
   */
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  /*
   void find_desired_criterion(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  bool always_positive() const;
  
protected:
  
  /**
   * @brief Copies the state of another PositiveSMCCriterion into this object.
   *
   * @param another The PositiveSMCCriterion instance to copy from.
   */
  void make_copy(const PositiveSMCCriterion &another);
  
};
}

#endif
