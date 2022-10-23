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
            SMCCriterion* criterion_in = NULL,
            SMCTermination* termination_in = NULL);
  
  //Sequencer(const std::vector< std::vector<double> > &schedules_in,
            //const std::vector<std::string> &variable_names_in);
  
  Sequencer(SMCWorker* the_worker_in,
            const std::vector< std::vector<double> > &schedules_in,
            const std::vector<std::string> &variable_names_in,
            SMCCriterion* criterion_in = NULL,
            SMCTermination* termination_in = NULL);

  Sequencer(const Sequencer &another);
  void operator=(const Sequencer &another);
  
  void find_desired_criterion(SMCOutput* current_state);
  
  arma::colvec find_next_target_bisection(SMCOutput* current_state,
                                          const Index* index);
  
  arma::colvec find_next_target_quantile(SMCOutput* current_state);
  
  void find_desired_criterion(SMCOutput* current_state,
                              const Parameters &conditioned_on_parameters);
  
  arma::colvec find_next_target_bisection(SMCOutput* current_state,
                                          const Index* index,
                                          const Parameters &conditioned_on_parameters);
  
  arma::colvec find_next_target_quantile(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters);
  
  arma::colvec subsample_find_next_target_bisection(SMCOutput* current_state,
                                                    const Index* index,
                                                    const Parameters &conditioned_on_parameters);
  
  arma::colvec subsample_find_next_target_quantile(SMCOutput* current_state,
                                                   const Parameters &conditioned_on_parameters);
  
  //void find_next_target_quantile(SMCOutput* current_state, const Parameters &conditioned_on_parameters);
  
  bool check_termination();
  
  void set_next_with_parameter(const Parameters &parameters_in);

protected:
  
  void setup(SMCWorker* the_worker_in,
             const std::vector< std::vector<double> > &schedules_in,
             const std::vector<std::string> &variable_names_in,
             SMCCriterion* criterion_in,
             SMCTermination* termination_in);

  void make_copy(const Sequencer &another);
  
  std::vector<double> current_values;
  
  // Creates a particular way of indexing some distributions.
  Mileometer mileometer;
  
  // The back end of the vector contains the variable we will cycle over fastest.
  std::vector< std::vector<double> > schedules;
  std::vector<std::string> variable_names;
  double direction; // does not need to be read in
  
  double current_score;
  
  double extra_bit;
  
  // not stored here
  SMCWorker* the_worker;
  
  // stored here
  SMCCriterion* criterion;
  SMCTermination* termination;

};

#endif
