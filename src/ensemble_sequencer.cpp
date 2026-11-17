#include "ensemble_sequencer.h"
#include "smc_criterion.h"
#include "smc_termination.h"
#include "ensemble_kalman_worker.h"
#include "ensemble_kalman_output.h"
#include "ensemble_kalman.h"

namespace ilike
{
EnsembleSequencer::EnsembleSequencer()
{
  this->criterion.reset();
  this->termination.reset();
  this->find_scale = false;
}

EnsembleSequencer::~EnsembleSequencer()
{
  this->criterion.reset();
  this->termination.reset();
}

EnsembleSequencer::EnsembleSequencer(EnsembleKalmanWorker* the_worker_in,
                                     const std::vector<double> &schedule_in,
                                     const std::string &variable_in,
                                     size_t number_of_bisections_in,
                                     SMCCriterion* criterion_in,
                                     SMCTermination* termination_in)
{
  this->setup(the_worker_in,
              schedule_in,
              variable_in,
              number_of_bisections_in,
              criterion_in,
              termination_in);
}

void EnsembleSequencer::setup(EnsembleKalmanWorker* the_worker_in,
                              const std::vector<double> &schedule_in,
                              const std::string &variable_in,
                              size_t number_of_bisections_in,
                              SMCCriterion* criterion_in,
                              SMCTermination* termination_in)
{
  this->the_worker = the_worker_in;
  
  this->criterion.reset(criterion_in);
  this->termination.reset(termination_in);
  this->variable_name = variable_in;
  this->number_of_bisections = number_of_bisections_in;
  
  this->set_schedule(schedule_in);
  
  this->find_scale = false;
}

void EnsembleSequencer::set_schedule(const std::vector<double> &schedule_in)
{
  this->schedule = schedule_in;
  
  if (this->schedule.size()<2)
    Rcpp::stop("EnsembleSequencer::setup - invalid schedule.");
  
  std::vector<size_t> sizes;
  sizes.reserve(1);
  sizes.push_back(this->schedule.size());
  this->use_final = true;
  
  this->mileometer = Mileometer(sizes);
  this->schedule_difference = arma::datum::inf;
  
  this->reset();
}

void EnsembleSequencer::set_initial_schedule_parameters()
{
  this->mileometer.reset();
  this->previous_bisect_value = this->schedule.front();
  this->mileometer.increment();
  
  double value_to_use;
  if ((this->use_final==true) && (this->mileometer[0]==0))
  {
    value_to_use = this->schedule[this->schedule.size()-1];
  }
  else
  {
    value_to_use = this->schedule[this->mileometer[0]];
  }
  this->schedule_parameters[this->variable_name] = value_to_use;
  this->current_bisect_value = value_to_use;
}

void EnsembleSequencer::set_schedule_parameters()
{
  double value_to_use;
  if ((this->use_final==true) && (this->mileometer[0]==0))
  {
    value_to_use = this->schedule[this->schedule.size()-1];
  }
  else
  {
    value_to_use = this->schedule[this->mileometer[0]];
  }
  this->schedule_parameters[this->variable_name] = value_to_use;
  
  this->mileometer.increment();
  
  this->current_bisect_value = this->schedule[this->mileometer[0]];
}

EnsembleSequencer::EnsembleSequencer(const EnsembleSequencer &another)
{
  this->make_copy(another);
}

EnsembleSequencer& EnsembleSequencer::operator=(const EnsembleSequencer &another)
{
  if(this == &another)
    return *this;
  
  this->criterion.reset();
  this->termination.reset();
  
  this->make_copy(another);
  
  return *this;
}

void EnsembleSequencer::make_copy(const EnsembleSequencer &another)
{
  this->schedule = another.schedule;
  this->direction = another.direction;
  this->current_bisect_value = another.current_bisect_value;
  this->previous_bisect_value = another.previous_bisect_value;
  this->mileometer = another.mileometer;
  this->current_score = another.current_score;
  this->number_of_bisections = another.number_of_bisections;
  this->the_worker = another.the_worker;
  if (another.criterion!=NULL)
    this->criterion.reset(another.criterion->duplicate());
  else
    this->criterion.reset();
  if (another.termination!=NULL)
    this->termination.reset(another.termination->duplicate());
  else
    this->termination.reset();
  this->variable_name = another.variable_name;
  this->schedule_difference = another.schedule_difference;
  this->schedule_parameters = another.schedule_parameters;
  this->use_final = another.use_final;
  this->find_scale = another.find_scale;
  this->scale_variable = another.scale_variable;
}

EnsembleSequencer::EnsembleSequencer(EnsembleSequencer &&another)
{
  this->make_copy(std::move(another));
}

EnsembleSequencer& EnsembleSequencer::operator=(EnsembleSequencer &&another)
{
  if(this == &another)
    return *this;
  
  this->criterion.reset();
  this->termination.reset();
  
  this->make_copy(std::move(another));
  
  return *this;
}

void EnsembleSequencer::make_copy(EnsembleSequencer &&another)
{
  this->schedule = std::move(another.schedule);
  this->direction = std::move(another.direction);
  this->current_bisect_value = std::move(another.current_bisect_value);
  this->previous_bisect_value = std::move(another.previous_bisect_value);
  this->mileometer = std::move(another.mileometer);
  this->current_score = std::move(another.current_score);
  this->number_of_bisections = std::move(another.number_of_bisections);
  this->the_worker = std::move(another.the_worker);
  this->criterion = std::move(another.criterion);
  this->termination = std::move(another.termination);
  this->variable_name = std::move(another.variable_name);
  this->schedule_difference = std::move(another.schedule_difference);
  this->schedule_parameters = std::move(another.schedule_parameters);
  this->use_final = std::move(another.use_final);
  this->find_scale = std::move(another.find_scale);
  this->scale_variable = std::move(another.scale_variable);
  
  another.schedule = std::vector<double>();
  another.direction = 0.0;
  another.current_bisect_value = 0.0;
  another.previous_bisect_value = 0.0;
  another.mileometer = Mileometer();
  another.current_score = 0.0;
  another.number_of_bisections = 0;
  another.the_worker = NULL;
  another.variable_name = "";
  another.schedule_difference = 0.0;
  another.schedule_parameters = Parameters();
  another.use_final = false;
  another.find_scale = false;
  another.scale_variable = "";
}

void EnsembleSequencer::find_desired_criterion(EnsembleKalmanOutput* current_state)
{
  this->criterion->find_desired_criterion(current_state);
}

// new method:
// sets desired ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

void EnsembleSequencer::find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                   const Index* index)
{
  if (current_state->skip_to_end_of_sequence)
  {
    if (this->current_bisect_value!=this->schedule.back())
      this->schedule_difference = this->schedule.back()-this->current_bisect_value;
    else
      this->schedule_difference = 1.0;
    
    this->current_bisect_value = this->schedule.back();
    this->schedule_parameters[this->variable_name] = this->current_bisect_value;
    this->mileometer.reset();
    
    return;
  }
  
  // used here since if we are not doing things adaptively, we never need to calculate a weight
  if (this->criterion->always_positive())
  {
    this->schedule_difference = this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1];
    this->current_bisect_value = this->schedule[this->mileometer.back()];
    this->set_schedule_parameters();
    
    return;
  }
  
  this->schedule_parameters[this->variable_name] = this->current_bisect_value;
  
  double small_value = std::min(abs((this->schedule[this->mileometer.back()] - this->schedule[this->mileometer.back()-1])),1.0)*pow(0.5,30.0);
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->direction * this->current_bisect_value>= this->direction * this->schedule[this->mileometer.back()] - small_value)
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(this->schedule[this->mileometer.back()]-this->previous_bisect_value));
    
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion.get())(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      this->schedule_difference = this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1];

      this->set_schedule_parameters();

      return;
    }
    
    this->current_bisect_value = this->previous_bisect_value;
  }

  this->previous_bisect_value = this->current_bisect_value;
  
  double next_value = this->schedule[this->mileometer.back()];
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  double starting_value = this->current_bisect_value;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_name] = new_bisect_value;
    
    if (new_bisect_value<starting_value)
    {
      this->current_bisect_value = starting_value;
      break;
    }
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-starting_value));
    
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion.get())(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }

  this->schedule_difference = this->current_bisect_value - starting_value;
  this->previous_bisect_value = this->current_bisect_value;
  
  if ((this->direction * this->current_bisect_value >= this->direction * this->schedule[this->mileometer.back()]-small_value))
  {
    this->mileometer.increment();
  }
  
}

arma::colvec EnsembleSequencer::find_next_target_quantile(EnsembleKalmanOutput* current_state)
{
  return arma::colvec();
}

void EnsembleSequencer::subsample_find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                             const Index* index)
{
  if (current_state->skip_to_end_of_sequence)
  {
    if (this->current_bisect_value!=this->schedule.back())
      this->schedule_difference = this->schedule.back()-this->schedule[this->mileometer.back()-1];
    else
      this->schedule_difference = 0.0;
    
    this->current_bisect_value = this->schedule.back();
    current_state->back().set_temperature(this->current_bisect_value);
    this->schedule_parameters[this->variable_name] = this->current_bisect_value;
    this->mileometer.reset();
    
    return;
  }
  
  // used here since if we are not doing things adaptively, we never need to calculate a weight
  if (this->criterion->always_positive())
  {
    
    this->schedule_difference = this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1];
    this->current_bisect_value = this->schedule[this->mileometer.back()];
    this->set_schedule_parameters();
    
    return;
  }
  
  this->schedule_parameters[this->variable_name] = this->current_bisect_value;
  
  double small_value = std::min(abs((this->schedule[this->mileometer.back()] - this->schedule[this->mileometer.back()-1])),1.0)*pow(0.5,30.0);
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->direction * this->current_bisect_value>= this->direction * this->schedule[this->mileometer.back()] - small_value)
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(this->schedule[this->mileometer.back()]-this->previous_bisect_value));
    
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion.get())(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      this->schedule_difference = this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1];

      this->set_schedule_parameters();

      return;
    }
    
    this->current_bisect_value = this->previous_bisect_value;
  }

  this->previous_bisect_value = this->current_bisect_value;
  
  double next_value = this->schedule[this->mileometer.back()];
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  double starting_value = this->current_bisect_value;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_name] = new_bisect_value;
    
    if (new_bisect_value<starting_value)
      new_bisect_value = starting_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(new_bisect_value-starting_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion.get())(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }

  this->schedule_difference = this->current_bisect_value - starting_value;
  this->previous_bisect_value = this->current_bisect_value;
  
  if ((this->direction * this->current_bisect_value >= this->direction * this->schedule[this->mileometer.back()]-small_value))
  {
    this->mileometer.increment();
  }
  
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
  new_schedule.push_back(this->schedule[this->mileometer.back()]);
  new_schedule.push_back(parameters_in[this->variable_name][0]);
  this->schedule = new_schedule;
  this->mileometer.reset_final_dimension(2);
  if (this->schedule.back()>=this->schedule.front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }
  
  this->mileometer.increment();
}

void EnsembleSequencer::reset()
{
  if (this->schedule.back()>=this->schedule.front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }
  
  this->set_initial_schedule_parameters();
  
  // Use epsilon_doubling in abc.r if we need help finding a maximum value to begin with.
}
}
