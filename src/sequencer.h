#ifndef SEQUENCER_H
#define SEQUENCER_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <string>
#include <vector>
#include "particles.h"
#include "mileometer.h"

namespace ilike
{

  class SMCCriterion;
  class SMCTermination;
  class SMCWorker;
  class SMCOutput;
  class Index;
  class ImportanceSampler;
  class SMCMCMCMove;

  class Sequencer
  {
  public:
    Sequencer();
    virtual ~Sequencer();

    Sequencer(SMCWorker *the_worker_in,
              const std::vector<double> &schedule_in,
              const std::string &variable_in,
              size_t number_of_bisections_in,
              SMCCriterion *criterion_in = NULL,
              SMCTermination *termination_in = NULL);

    Sequencer(SMCWorker *the_worker_in,
              const std::vector<std::vector<double>> &schedules_in,
              const std::vector<std::string> &variable_names_in,
              size_t number_of_bisections_in,
              SMCCriterion *criterion_in = NULL,
              SMCTermination *termination_in = NULL);

    Sequencer(const Sequencer &another);
    Sequencer &operator=(const Sequencer &another);

    Sequencer(Sequencer &&another);
    Sequencer &operator=(Sequencer &&another);

    void find_desired_criterion(SMCOutput *current_state);

    // The worker this sequencer weights through (after its SMC has been
    // copied, B29).
    void set_worker(SMCWorker *the_worker_in) { this->the_worker = the_worker_in; }

    void find_next_target_bisection(SMCOutput *current_state,
                                    const Index *index);

    void find_next_target_quantile(SMCOutput *current_state,
                                   const Index *index);

    void subsample_find_desired_criterion(SMCOutput *current_state);

    void subsample_find_next_target_bisection(SMCOutput *current_state,
                                              const Index *index);

    void subsample_find_next_target_quantile(SMCOutput *current_state,
                                             const Index *index);

    bool check_termination();

    // Progress along the (last) schedule for a progress bar (dev/plan.md
    // 3.11): fraction in [0,1] of the way from its start to its end, or a
    // negative value when an end is infinite; status: "variable=value".
    double progress_fraction(std::string &status) const;

    void set_next_with_parameter(const Parameters &parameters_in);

    Parameters schedule_parameters;

    void reset();

    // (dev/plan.md 3.22 U4) parameters moved along the last schedule variable
    // t in [0, 1], from start to end (an annealed path between two points)
    void set_path(const std::vector<std::string> &names, const Parameters &start, const Parameters &end);
    // the last schedule variable's current value, and its first point; and
    // setting it (with a path's parameters following)
    double last_value() const { return this->schedule_parameters[this->variable_names.back()][0]; }
    double first_point() const { return this->schedules.back().front(); }
    void set_last_value(double value) { this->set_schedule_value(this->variable_names.back(), value); }

  protected:
    void set_schedule_value(const std::string &name, double value);
    void apply_path();
    std::vector<std::string> path_names;
    Parameters path_start;
    Parameters path_end;

    friend ImportanceSampler;
    friend SMCMCMCMove; // (dev/plan.md 3.21 (c4): the schedule, for early stopping)

    void setup(SMCWorker *the_worker_in,
               const std::vector<std::vector<double>> &schedules_in,
               const std::vector<std::string> &variable_names_in,
               size_t number_of_bisections_in,
               SMCCriterion *criterion_in,
               SMCTermination *termination_in);

    void set_initial_schedule_parameters();
    void set_schedule_parameters();

    void make_copy(const Sequencer &another);
    void make_copy(Sequencer &&another);

    // Creates a particular way of indexing some distributions.
    Mileometer mileometer;

    double current_bisect_value;
    double previous_bisect_value;

    // The back end of the vector contains the variable we will cycle over fastest.
    std::vector<std::vector<double>> schedules;
    std::vector<std::string> variable_names;
    std::vector<bool> use_final; // for each variable, if we reach the last value, do we set the parameters to use the final value, or the one at the beginning.
    double direction;            // does not need to be read in

    double current_score;

    double schedule_difference;

    size_t number_of_bisections;

    // info required to automatically determine target
    // so far, for EnKI-ABC only
    std::string scale_variable;
    bool find_scale;
    std::vector<std::string> scale_variables;

    // not stored here
    SMCWorker *the_worker;

    // stored here
    // owned (dev/plan.md 3.5, CH3)
    std::unique_ptr<SMCCriterion> criterion;
    std::unique_ptr<SMCTermination> termination;
  };
}

#endif
