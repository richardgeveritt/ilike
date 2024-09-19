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
  
  EnsembleSequencer();
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
  
  EnsembleSequencer(const EnsembleSequencer &another);
  EnsembleSequencer& operator=(const EnsembleSequencer &another);
  
  EnsembleSequencer(EnsembleSequencer &&another);
  EnsembleSequencer& operator=(EnsembleSequencer &&another);
  
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  void find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                  const Index* index);
  
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
  
  void set_schedule(const std::vector<double> &schedule_in);
  
  bool check_termination();
  
  void set_next_with_parameter(const Parameters &parameters_in);
  
  Parameters schedule_parameters;
  
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
  
  void set_initial_schedule_parameters();
  void set_schedule_parameters();
  
  void make_copy(const EnsembleSequencer &another);
  void make_copy(EnsembleSequencer &&another);
  
  //double current_value;
  
  // Creates a particular way of indexing some distributions.
  Mileometer mileometer;
  
  double current_bisect_value;
  double previous_bisect_value;
  
  // The back end of the vector contains the variable we will cycle over fastest.
  std::vector<double> schedule;
  std::string variable_name;
  bool use_final;
  double direction; // does not need to be read in
  
  double current_score;
  
  double schedule_difference;
  
  size_t number_of_bisections;
  
  
  // info required to automatically determine target
  // so far, for EnKI-ABC only
  std::string scale_variable;
  bool find_scale;
  //arma::colvec simulation_scale;
  
  //double extra_bit;
  
  // not stored here
  EnsembleKalmanWorker* the_worker;
  
  // stored here
  SMCCriterion* criterion;
  SMCTermination* termination;
  
};
}

#endif
