#include "ensemble_sequencer.h"
#include "smc_criterion.h"
#include "smc_termination.h"
#include "ensemble_kalman_worker.h"
#include "ensemble_kalman_output.h"
#include "ensemble_kalman.h"

EnsembleSequencer::EnsembleSequencer()
{
  this->criterion = NULL;
  this->termination = NULL;
}

EnsembleSequencer::~EnsembleSequencer()
{
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;
}

//EnsembleSequencer::EnsembleSequencer(const std::vector<double> &schedule_in,
//                     const std::string &variable_in)
//{
//  std::vector< std::vector<double> > schedules_in;
//  schedules_in.push_back(schedule_in);
//  std::vector<std::string> variable_names_in;
//  variable_names_in.push_back(variable_in);
//  SMCCriterion* criterion_in = NULL;
//  double desired_criterion_in = 0.0;
//  this->setup(schedules_in, variable_names_in, desired_criterion_in, criterion_in);
//}

EnsembleSequencer::EnsembleSequencer(EnsembleKalmanWorker* the_worker_in,
                                     const std::vector<double> &schedule_in,
                                     SMCCriterion* criterion_in,
                                     SMCTermination* termination_in)
{
  this->setup(the_worker_in,
              schedule_in,
              criterion_in,
              termination_in);
}

void EnsembleSequencer::setup(EnsembleKalmanWorker* the_worker_in,
                              const std::vector<double> &schedule_in,
                              SMCCriterion* criterion_in,
                              SMCTermination* termination_in)
{
  this->the_worker = the_worker_in;
  this->schedule = schedule_in;
  this->criterion = criterion_in;
  this->termination = termination_in;
  
  if (this->schedule.size()<2)
      throw std::runtime_error("EnsembleSequencer::setup - invalid schedule.");
      
  if (this->schedule.back()>=this->schedule.front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }
  
  this->extra_bit = 0.001*(this->schedule.back()-this->schedule.front());
  
  std::vector<size_t> sizes;
  sizes.reserve(1);
  this->current_value = this->schedule.front();
  sizes.push_back(this->schedule.size());
  
  this->mileometer = Mileometer(sizes);
  this->mileometer.increment();
  
  // Use epsilon_doubling in abc.r if we need help finding a maximum value to begin with.
}

EnsembleSequencer::EnsembleSequencer(const EnsembleSequencer &another)
{
  this->make_copy(another);
}

void EnsembleSequencer::operator=(const EnsembleSequencer &another)
{
  if(this == &another)
    return;
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;

  this->make_copy(another);
}

void EnsembleSequencer::make_copy(const EnsembleSequencer &another)
{
  this->schedule = another.schedule;
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = another.direction;
  this->current_value = another.current_value;
  this->mileometer = another.mileometer;
  this->current_score = another.current_score;
  this->extra_bit = another.extra_bit;
  this->the_worker = another.the_worker;
  if (another.criterion!=NULL)
    this->criterion = another.criterion->duplicate();
  else
    this->criterion = NULL;
  if (another.termination!=NULL)
    this->termination = another.termination->duplicate();
  else
    this->termination = NULL;
  this->variable_name = another.variable_name;
}

void EnsembleSequencer::find_desired_criterion(EnsembleKalmanOutput* current_state)
{
  this->criterion->find_desired_criterion(current_state);
}

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

void EnsembleSequencer::find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                   const Index* index)
{
  double current_bisect_value = this->current_value;
  double previous_bisect_value = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-current_bisect_value));
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    //this->the_worker->weight(current_state->all_ensembles.back(),
    //                         new_bisect_value-previous_bisect_value);
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    //this->mileometer.increment();
  }
  
  this->current_value = new_bisect_value;
  current_state->back().set_temperature(this->current_value);
}

arma::colvec EnsembleSequencer::find_next_target_quantile(EnsembleKalmanOutput* current_state)
{
  return arma::colvec();
}

void EnsembleSequencer::find_desired_criterion(EnsembleKalmanOutput* current_state,
                                       const Parameters &conditioned_on_parameters)
{
  this->criterion->find_desired_criterion(current_state, conditioned_on_parameters);
}

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

void EnsembleSequencer::find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                   const Index* index,
                                                   const Parameters &conditioned_on_parameters)
{
  double current_bisect_value = this->current_value;
  double previous_bisect_value = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-current_bisect_value),
                                                  conditioned_on_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    //this->the_worker->weight(current_state->all_ensembles.back(),
    //                         new_bisect_value-previous_bisect_value);
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    //this->mileometer.increment();
  }
  
  this->current_value = new_bisect_value;
  current_state->back().set_temperature(this->current_value);
}

arma::colvec EnsembleSequencer::find_next_target_quantile(EnsembleKalmanOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  return arma::colvec();
}

void EnsembleSequencer::subsample_find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                             const Index* index,
                                                             const Parameters &conditioned_on_parameters)
{
  double current_bisect_value = this->current_value;
  double previous_bisect_value = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(new_bisect_value-current_bisect_value),
                                                            conditioned_on_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    //this->the_worker->weight(current_state->all_ensembles.back(),
    //                         new_bisect_value-previous_bisect_value);
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    //this->mileometer.increment();
  }
  
  this->current_value = new_bisect_value;
  current_state->back().set_temperature(this->current_value);
}

arma::colvec EnsembleSequencer::subsample_find_next_target_quantile(EnsembleKalmanOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  return arma::colvec();
}

bool EnsembleSequencer::check_termination()
{
  if (this->mileometer.at_start())
    return true;
  else
  {
    if ( (this->termination!=NULL) && ( this->termination->terminate(this->current_score) ) )
      return true;
  }
  return false;
}

void EnsembleSequencer::set_next_with_parameter(const Parameters &parameters_in)
{
  std::vector<double> new_schedule;
  new_schedule.reserve(2);
  new_schedule.push_back(this->current_value);
  new_schedule.push_back(parameters_in[this->variable_name][0]);
  this->schedule = new_schedule;
  this->mileometer.current_index.back() = 0;
  this->mileometer.limits.back() = 2;
  
  this->mileometer.increment();
}
