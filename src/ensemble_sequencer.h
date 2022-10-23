#ifndef ENSEMBLESEQUENCER_H
#define ENSEMBLESEQUENCER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "particles.h"
#include "mileometer.h"

class SMCCriterion;
class SMCTermination;
class EnsembleKalmanWorker;
class EnsembleKalmanOutput;
class IterativeEnsembleKalmanInversion;
class EnsembleKalmanMFDS;

class EnsembleSequencer
{
public:

  EnsembleSequencer();
  virtual ~EnsembleSequencer();
  
  //EnsembleSequencer(const std::vector<double> &schedule_in,
            //const std::string &variable_in);
  
  EnsembleSequencer(EnsembleKalmanWorker* the_worker_in,
                    const std::vector<double> &schedule_in,
                    SMCCriterion* criterion_in,
                    SMCTermination* termination_in);
  
  //EnsembleSequencer(const std::vector< std::vector<double> > &schedules_in,
            //const std::vector<std::string> &variable_names_in);

  EnsembleSequencer(const EnsembleSequencer &another);
  void operator=(const EnsembleSequencer &another);
  
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  void find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                  const Index* index);
  
  arma::colvec find_next_target_quantile(EnsembleKalmanOutput* current_state);
  
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
  
  //void find_next_target_quantile(SMCOutput* current_state, const Parameters &conditioned_on_parameters);
  
  bool check_termination();
  
  void set_next_with_parameter(const Parameters &parameters_in);

protected:
  
  friend IterativeEnsembleKalmanInversion;
  friend EnsembleKalmanMFDS;
  
  void setup(EnsembleKalmanWorker* the_worker_in,
             const std::vector<double> &schedules_in,
             SMCCriterion* criterion_in,
             SMCTermination* termination_in);

  void make_copy(const EnsembleSequencer &another);
  
  double current_value;
  
  // Creates a particular way of indexing some distributions.
  Mileometer mileometer;
  
  // The back end of the vector contains the variable we will cycle over fastest.
  std::vector<double> schedule;
  std::string variable_name;
  double direction; // does not need to be read in
  
  double current_score;
  
  double extra_bit;
  
  // not stored here
  EnsembleKalmanWorker* the_worker;
  
  // stored here
  SMCCriterion* criterion;
  SMCTermination* termination;

};

#endif
