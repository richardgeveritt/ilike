#ifndef SEQUENCER_H
#define SEQUENCER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "particles.h"
#include "mileometer.h"

class SMCCriterion;
class SMCTermination;
class SMCWorker;
class SMCOutput;
class Index;

class Sequencer
{
public:

  Sequencer();
  virtual ~Sequencer();
  
  //Sequencer(const std::vector<double> &schedule_in,
            //const std::string &variable_in);
  
  Sequencer(SMCWorker* the_worker_in,
            const std::vector<double> &schedule_in,
            const std::string &variable_in,
            size_t number_of_bisections_in,
            SMCCriterion* criterion_in = NULL,
            SMCTermination* termination_in = NULL);
  
  //Sequencer(const std::vector< std::vector<double> > &schedules_in,
            //const std::vector<std::string> &variable_names_in);
  
  Sequencer(SMCWorker* the_worker_in,
            const std::vector< std::vector<double> > &schedules_in,
            const std::vector<std::string> &variable_names_in,
            size_t number_of_bisections_in,
            SMCCriterion* criterion_in = NULL,
            SMCTermination* termination_in = NULL);

  Sequencer(const Sequencer &another);
  Sequencer& operator=(const Sequencer &another);
  
  Sequencer(Sequencer &&another);
  Sequencer& operator=(Sequencer &&another);
  
  void find_desired_criterion(SMCOutput* current_state);
  
  void find_next_target_bisection(SMCOutput* current_state,
                                  const Index* index);
  
  void find_next_target_quantile(SMCOutput* current_state,
                                 const Index* index);
  
  void subsample_find_desired_criterion(SMCOutput* current_state);
  
  void subsample_find_next_target_bisection(SMCOutput* current_state,
                                            const Index* index);
  
  void subsample_find_next_target_quantile(SMCOutput* current_state,
                                           const Index* index);
  
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
  
  //void find_next_target_quantile(SMCOutput* current_state, const Parameters &conditioned_on_parameters);
  
  bool check_termination();
  
  void set_next_with_parameter(const Parameters &parameters_in);
  
  Parameters schedule_parameters;
  
  void reset();

protected:
  
  void setup(SMCWorker* the_worker_in,
             const std::vector< std::vector<double> > &schedules_in,
             const std::vector<std::string> &variable_names_in,
             size_t number_of_bisections_in,
             SMCCriterion* criterion_in,
             SMCTermination* termination_in);
  
  void set_initial_schedule_parameters();
  void set_schedule_parameters();

  void make_copy(const Sequencer &another);
  void make_copy(Sequencer &&another);
  
  //std::vector<double> current_values;
  
  // Creates a particular way of indexing some distributions.
  Mileometer mileometer;
  
  double current_bisect_value;
  double previous_bisect_value;
  
  // The back end of the vector contains the variable we will cycle over fastest.
  std::vector< std::vector<double> > schedules;
  std::vector<std::string> variable_names;
  std::vector<bool> use_final; // for each variable, if we reach the last value, do we set the parameters to use the final value, or the one at the beginning.
  double direction; // does not need to be read in
  
  double current_score;
  
  double schedule_difference;
  
  size_t number_of_bisections;
  
  //double extra_bit;
  
  // not stored here
  SMCWorker* the_worker;
  
  // stored here
  SMCCriterion* criterion;
  SMCTermination* termination;

};

#endif
