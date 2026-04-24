#ifndef SEQUENCER_H
#define SEQUENCER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "particles.h"
#include "mileometer.h"

namespace ilike
{
  /**
   * @file sequencer.h
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
  class SMCWorker;
  class SMCOutput;
  class Index;
  class ImportanceSampler;

  class Sequencer
  {
  public:
    /**
     * @brief Performs the sequencer operation.
     */
    Sequencer();
    /**
     * @brief Performs the ~sequencer operation.
     */
    virtual ~Sequencer();

    // Sequencer(const std::vector<double> &schedule_in,
    // const std::string &variable_in);

    Sequencer(SMCWorker *the_worker_in,
              const std::vector<double> &schedule_in,
              const std::string &variable_in,
              size_t number_of_bisections_in,
              SMCCriterion *criterion_in = NULL,
              SMCTermination *termination_in = NULL);

    // Sequencer(const std::vector< std::vector<double> > &schedules_in,
    // const std::vector<std::string> &variable_names_in);

    Sequencer(SMCWorker *the_worker_in,
              const std::vector<std::vector<double>> &schedules_in,
              const std::vector<std::string> &variable_names_in,
              size_t number_of_bisections_in,
              SMCCriterion *criterion_in = NULL,
              SMCTermination *termination_in = NULL);

    /**
     * @brief Performs the sequencer operation.
     *
     * @param another The SMCCriterion instance to copy from.
     */
    Sequencer(const Sequencer &another);
    Sequencer &operator=(const Sequencer &another);

    /**
     * @brief Performs the sequencer operation.
     *
     * @param another The SMCCriterion instance to copy from.
     */
    Sequencer(Sequencer &&another);
    Sequencer &operator=(Sequencer &&another);

    /**
     * @brief Finds desired criterion.
     *
     * @param current_state The current state.
     */
    void find_desired_criterion(SMCOutput *current_state);

    void find_next_target_bisection(SMCOutput *current_state,
                                    const Index *index);

    void find_next_target_quantile(SMCOutput *current_state,
                                   const Index *index);

    /**
     * @brief Performs the subsample find desired criterion operation.
     *
     * @param current_state The current state.
     */
    void subsample_find_desired_criterion(SMCOutput *current_state);

    void subsample_find_next_target_bisection(SMCOutput *current_state,
                                              const Index *index);

    void subsample_find_next_target_quantile(SMCOutput *current_state,
                                             const Index *index);

    /*
     void find_desired_criterion(SMCOutput* current_state,
     const Parameters &conditioned_on_parameters);

     void subsample_find_desired_criterion(SMCOutput* current_state,
     const Parameters &conditioned_on_parameters);
     */

    /*
     void find_next_target_bisection(SMCOutput* current_state,
     const Index* index,
     const Parameters &conditioned_on_parameters);

     void find_next_target_quantile(SMCOutput* current_state,
     const Index* index,
     const Parameters &conditioned_on_parameters);

     void subsample_find_next_target_bisection(SMCOutput* current_state,
     const Index* index,
     const Parameters &conditioned_on_parameters);

     void subsample_find_next_target_quantile(SMCOutput* current_state,
     const Index* index,
     const Parameters &conditioned_on_parameters);
     */

    // void find_next_target_quantile(SMCOutput* current_state, const Parameters &conditioned_on_parameters);

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
    friend ImportanceSampler;

    void setup(SMCWorker *the_worker_in,
               const std::vector<std::vector<double>> &schedules_in,
               const std::vector<std::string> &variable_names_in,
               size_t number_of_bisections_in,
               SMCCriterion *criterion_in,
               SMCTermination *termination_in);

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
    void make_copy(const Sequencer &another);
    /**
     * @brief Copies the state of another SMCCriterion into this object.
     *
     * @param another The SMCCriterion instance to copy from.
     */
    void make_copy(Sequencer &&another);

    // std::vector<double> current_values;

    // Creates a particular way of indexing some distributions.
    /** @brief The mileometer. */
    Mileometer mileometer;

    /** @brief The current bisect value. */
    double current_bisect_value;
    /** @brief The previous bisect value. */
    double previous_bisect_value;

    // The back end of the vector contains the variable we will cycle over fastest.
    /** @brief The schedules. */
    std::vector<std::vector<double>> schedules;
    /** @brief The variable names. */
    std::vector<std::string> variable_names;
    /** @brief The use final. */
    std::vector<bool> use_final; // for each variable, if we reach the last value, do we set the parameters to use the final value, or the one at the beginning.
    /** @brief The direction. */
    double direction;            // does not need to be read in

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
    /** @brief The scale variables. */
    std::vector<std::string> scale_variables;

    // double extra_bit;

    // not stored here
    SMCWorker *the_worker;

    // stored here
    SMCCriterion *criterion;
    SMCTermination *termination;
  };
}

#endif
