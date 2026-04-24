#ifndef ENSEMBLESEQUENCER_H
#define ENSEMBLESEQUENCER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "particles.h"
#include "mileometer.h"

namespace ilike
{
  /**
   * @file ensemble_sequencer.h
   * @brief Defines the SMCCriterion class.
   *
   * Implements an SMC resampling/adaptation criterion based on default. Used by the SMC algorithm to decide when to resample particles.
   *
   * @namespace ilike
   * @class SMCCriterion
   * @brief The smc criterion class.
   */


class SMCCriterion;
class SMCTermination;
class EnsembleKalmanWorker;
class EnsembleKalmanOutput;
class EnsembleKalmanInversion;
class EnsembleKalmanMFDS;
class EnsembleKalmanFilter;
class EnsembleKalman;

class EnsembleSequencer
{
public:
  
  /**
   * @brief Performs the ensemblesequencer operation.
   */
  EnsembleSequencer();
  /**
   * @brief Performs the ~ensemblesequencer operation.
   */
  virtual ~EnsembleSequencer();
  
  //EnsembleSequencer(const std::vector<double> &schedule_in,
  //const std::string &variable_in);
  
  EnsembleSequencer(EnsembleKalmanWorker* the_worker_in,
                    const std::vector<double> &schedule_in,
                    const std::string &variable_in,
                    size_t number_of_bisections_in,
                    SMCCriterion* criterion_in = NULL,
                    SMCTermination* termination_in = NULL);
  
  //EnsembleSequencer(const std::vector< std::vector<double> > &schedules_in,
  //const std::vector<std::string> &variable_names_in);
  
  /**
   * @brief Performs the ensemblesequencer operation.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  EnsembleSequencer(const EnsembleSequencer &another);
  /**
   * @brief Assignment operator for SMCCriterion.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  EnsembleSequencer& operator=(const EnsembleSequencer &another);
  
  /**
   * @brief Performs the ensemblesequencer operation.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  EnsembleSequencer(EnsembleSequencer &&another);
  /**
   * @brief Assignment operator for SMCCriterion.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  EnsembleSequencer& operator=(EnsembleSequencer &&another);
  
  /**
   * @brief Finds desired criterion.
   *
   * @param current_state The current state.
   */
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  void find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                  const Index* index);
  
  /**
   * @brief Finds next target quantile.
   *
   * @param current_state The current state.
   *
   * @return The result.
   */
  arma::colvec find_next_target_quantile(EnsembleKalmanOutput* current_state);
  
  void subsample_find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                            const Index* index);
  
  /*
   void find_desired_criterion(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   
   void find_next_target_bisection(EnsembleKalmanOutput* current_state,
   const Index* index,
   const Parameters &conditioned_on_parameters);
   
   arma::colvec find_next_target_quantile(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   
   void subsample_find_next_target_bisection(EnsembleKalmanOutput* current_state,
   const Index* index,
   const Parameters &conditioned_on_parameters);
   
   arma::colvec subsample_find_next_target_quantile(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  //void find_next_target_quantile(SMCOutput* current_state, const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Sets the schedule.
   *
   * @param schedule_in The schedule.
   */
  void set_schedule(const std::vector<double> &schedule_in);
  
  /**
   * @brief Performs the check termination operation.
   *
   * @return The result.
   */
  bool check_termination();
  
  /**
   * @brief Sets the next with parameter.
   *
   * @param parameters_in The parameters.
   */
  void set_next_with_parameter(const Parameters &parameters_in);
  
  Parameters schedule_parameters;
  
  /**
   * @brief Performs the reset operation.
   */
  void reset();
  
protected:
  
  friend EnsembleKalmanInversion;
  friend EnsembleKalmanMFDS;
  friend EnsembleKalmanFilter;
  friend EnsembleKalman;
  
  void setup(EnsembleKalmanWorker* the_worker_in,
             const std::vector<double> &schedules_in,
             const std::string &variable_in,
             size_t number_of_bisections_in,
             SMCCriterion* criterion_in,
             SMCTermination* termination_in);
  
  /**
   * @brief Sets the initial schedule parameters.
   */
  void set_initial_schedule_parameters();
  /**
   * @brief Sets the schedule parameters.
   */
  void set_schedule_parameters();
  
  /**
   * @brief Copies the state of another SMCCriterion into this object.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  void make_copy(const EnsembleSequencer &another);
  /**
   * @brief Copies the state of another SMCCriterion into this object.
   *
   * @param another The SMCCriterion instance to copy from.
   */
  void make_copy(EnsembleSequencer &&another);
  
  //double current_value;
  
  // Creates a particular way of indexing some distributions.
  /** @brief The mileometer. */
  Mileometer mileometer;
  
  /** @brief The current bisect value. */
  double current_bisect_value;
  /** @brief The previous bisect value. */
  double previous_bisect_value;
  
  // The back end of the vector contains the variable we will cycle over fastest.
  /** @brief The schedule. */
  std::vector<double> schedule;
  /** @brief The variable name. */
  std::string variable_name;
  /** @brief The use final. */
  bool use_final;
  /** @brief The direction. */
  double direction; // does not need to be read in
  
  /** @brief The current score. */
  double current_score;
  
  /** @brief The schedule difference. */
  double schedule_difference;
  
  /** @brief The number of bisections. */
  size_t number_of_bisections;
  
  
  // info required to automatically determine target
  // so far, for EnKI-ABC only
  /** @brief The scale variable. */
  std::string scale_variable;
  /** @brief The find scale. */
  bool find_scale;
  //arma::colvec simulation_scale;
  
  //double extra_bit;
  
  // not stored here
  /** @brief The the worker. */
  EnsembleKalmanWorker* the_worker;
  
  // stored here
  /** @brief The criterion. */
  SMCCriterion* criterion;
  /** @brief The termination. */
  SMCTermination* termination;
  
};
}

#endif
