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
  this->schedule = schedule_in;
  this->criterion = criterion_in;
  this->termination = termination_in;
  this->variable_name = variable_in;
  this->number_of_bisections = number_of_bisections_in;
  
  if (this->schedule.size()<2)
      Rcpp::stop("EnsembleSequencer::setup - invalid schedule.");
  
  std::vector<size_t> sizes;
  sizes.reserve(1);
  //this->current_value = this->schedule.front();
  sizes.push_back(this->schedule.size());
  this->use_final = true;

  this->mileometer = Mileometer(sizes);
      
  this->reset();
}

void EnsembleSequencer::set_initial_schedule_parameters()
{
  this->mileometer.reset();
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
  
  this->current_bisect_value = value_to_use;
}

EnsembleSequencer::EnsembleSequencer(const EnsembleSequencer &another)
{
  this->make_copy(another);
}

EnsembleSequencer& EnsembleSequencer::operator=(const EnsembleSequencer &another)
{
  if(this == &another)
    return *this;
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;

  this->make_copy(another);
  
  return *this;
}

void EnsembleSequencer::make_copy(const EnsembleSequencer &another)
{
  this->schedule = another.schedule;
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = another.direction;
  this->current_bisect_value = another.current_bisect_value;
  this->mileometer = another.mileometer;
  this->current_score = another.current_score;
  this->number_of_bisections = another.number_of_bisections;
  //this->extra_bit = another.extra_bit;
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
  this->schedule_parameters = another.schedule_parameters;
  this->use_final = another.use_final;
}

EnsembleSequencer::EnsembleSequencer(EnsembleSequencer &&another)
{
  this->make_copy(std::move(another));
}

EnsembleSequencer& EnsembleSequencer::operator=(EnsembleSequencer &&another)
{
  if(this == &another)
    return *this;
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;
  
  this->make_copy(std::move(another));
  
  return *this;
}

void EnsembleSequencer::make_copy(EnsembleSequencer &&another)
{
  this->schedule = std::move(another.schedule);
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = std::move(another.direction);
  this->current_bisect_value = std::move(another.current_bisect_value);
  this->mileometer = std::move(another.mileometer);
  this->current_score = std::move(another.current_score);
  this->number_of_bisections = std::move(another.number_of_bisections);
  //this->extra_bit = another.extra_bit;
  this->the_worker = std::move(another.the_worker);
  if (another.criterion!=NULL)
    this->criterion = another.criterion;
  else
    this->criterion = NULL;
  if (another.termination!=NULL)
    this->termination = another.termination;
  else
    this->termination = NULL;
  this->variable_name = std::move(another.variable_name);
  this->schedule_parameters = std::move(another.schedule_parameters);
  this->use_final = std::move(another.use_final);
  
  another.schedule = std::vector<double>();
  another.direction = 0.0;
  another.current_bisect_value = 0.0;
  another.mileometer = Mileometer();
  another.current_score = 0.0;
  another.number_of_bisections = 0;
  another.the_worker = NULL;
  another.criterion = NULL;
  another.termination = NULL;
  another.variable_name = "";
  another.schedule_parameters = Parameters();
  another.use_final = false;
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
  if (this->criterion->always_positive())
  {
    //this->mileometer.increment();
    this->set_schedule_parameters();
    current_state->back().set_temperature(this->current_bisect_value);
    return;
  }
  
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->current_bisect_value>=this->schedule[this->mileometer.back()])
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1]));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      //this->mileometer.increment();
      this->set_schedule_parameters();
      current_state->back().set_temperature(this->current_bisect_value);
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
    
    this->current_bisect_value = this->schedule[this->mileometer.back()-1];
  }
  else
  {
    this->schedule_parameters[this->variable_name] = this->schedule[this->mileometer.back()];
    
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(this->schedule[this->mileometer.back()]-this->current_bisect_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      current_state->back().set_temperature(this->schedule[this->mileometer.back()]);
      //this->mileometer.increment();
      this->set_schedule_parameters();
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
  }
  
  double starting_bisect_value = this->current_bisect_value;
  double next_value = this->schedule[this->mileometer.back()];
  // target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_name] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-starting_bisect_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  current_state->back().set_temperature(this->current_bisect_value);
  
  /*
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_value;
  //double previous_bisect_value;// = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  //arma::colvec incremental_log_weights;
  
  //arma::colvec incremental_log_weights;
  
  if (this->criterion->always_positive())
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                index,
                                                1.0/(target_values.back()-current_bisect_value));
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->all_ensembles.back());
  
  if (this->current_score>=0.0)
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  double starting_bisect_value = current_bisect_value;
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-starting_bisect_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    //previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_value = new_bisect_value;
  this->schedule_parameters[this->variable_name] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    //this->the_worker->weight(current_state->all_ensembles.back(),
    //                         new_bisect_value-previous_bisect_value);
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
  }
  
  current_state->back().set_temperature(this->current_value);
  */
}

arma::colvec EnsembleSequencer::find_next_target_quantile(EnsembleKalmanOutput* current_state)
{
  return arma::colvec();
}

void EnsembleSequencer::subsample_find_next_target_bisection(EnsembleKalmanOutput* current_state,
                                                             const Index* index)
{
  if (this->criterion->always_positive())
  {
    //this->mileometer.increment();
    this->set_schedule_parameters();
    current_state->back().set_temperature(this->current_bisect_value);
    return;
  }
  
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->current_bisect_value>=this->schedule[this->mileometer.back()])
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                  index,
                                                  1.0/(this->schedule[this->mileometer.back()]-this->schedule[this->mileometer.back()-1]));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      //this->mileometer.increment();
      this->set_schedule_parameters();
      current_state->back().set_temperature(this->current_bisect_value);
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
    
    this->current_bisect_value = this->schedule[this->mileometer.back()-1];
  }
  else
  {
    this->schedule_parameters[this->variable_name] = this->schedule[this->mileometer.back()];
    
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(this->schedule[this->mileometer.back()]-this->current_bisect_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    if (this->current_score>=0.0)
    {
      current_state->back().set_temperature(this->schedule[this->mileometer.back()]);
      //this->mileometer.increment();
      this->set_schedule_parameters();
      //current_state->back().set_temperature(this->current_bisect_value);
      //current_state->back().set_temperature(this->schedule[this->mileometer.back()]);
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
  }
  
  double starting_bisect_value = this->current_bisect_value;
  double next_value = this->schedule[this->mileometer.back()];
  // target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_name] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-starting_bisect_value));
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  current_state->back().set_temperature(this->current_bisect_value);
  
  /*
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_value;
  //double previous_bisect_value;// = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  //arma::colvec incremental_log_weights;
  
  if (this->criterion->always_positive())
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  // thing to sort out:
  // don't actually want the last target value - want to reach the next one from the current
  // when we reach the target, want to increment, and set the current to be this
  
  this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                          index,
                                                          1.0/(target_values.back()-current_bisect_value));
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->all_ensembles.back());
  
  if (this->current_score>=0.0)
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  double starting_bisect_value = current_bisect_value;
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(new_bisect_value-starting_bisect_value));
    
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    //previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_value = new_bisect_value;
  this->schedule_parameters[this->variable_name] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
  }
  
  current_state->back().set_temperature(this->current_value);
  */
}

/*
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
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_value;
  //double previous_bisect_value;// = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  //arma::colvec incremental_log_weights;
  
  //arma::colvec incremental_log_weights;
  
  if (this->criterion->always_positive())
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                index,
                                                1.0/(target_values.back()-current_bisect_value),
                                                conditioned_on_parameters);
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->all_ensembles.back());
  
  if (this->current_score>=0.0)
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  double starting_bisect_value = current_bisect_value;
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->weight(&current_state->back(),
                                                  index,
                                                  1.0/(new_bisect_value-starting_bisect_value),
                                                  conditioned_on_parameters);
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    //previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_value = new_bisect_value;
  this->schedule_parameters[this->variable_name] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
  }
  
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
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_value;
  //double previous_bisect_value;// = this->schedule.front();
  std::vector< std::vector<double> > schedules;
  schedules.push_back(this->schedule);
  std::vector<double> target_values = this->mileometer.get_current_values(schedules);
  //arma::colvec incremental_log_weights;
  
  if (this->criterion->always_positive())
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  // thing to sort out:
  // don't actually want the last target value - want to reach the next one from the current
  // when we reach the target, want to increment, and set the current to be this
  
  this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                          index,
                                                          1.0/(target_values.back()-current_bisect_value),
                                                          conditioned_on_parameters);
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->all_ensembles.back());
  
  if (this->current_score>=0.0)
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
    current_state->back().set_temperature(this->current_value);
    
    return;
  }
  
  double starting_bisect_value = current_bisect_value;
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_enk->the_worker->subsample_weight(&current_state->back(),
                                                            index,
                                                            1.0/(new_bisect_value-starting_bisect_value),
                                                            conditioned_on_parameters);
    
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->all_ensembles.back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    //previous_bisect_value = current_bisect_value;
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_value = new_bisect_value;
  this->schedule_parameters[this->variable_name] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    this->schedule_parameters[this->variable_name] = target_values.back();
    
    this->mileometer.increment();
    
    this->current_value = target_values.back();
  }
  
  current_state->back().set_temperature(this->current_value);
}

arma::colvec EnsembleSequencer::subsample_find_next_target_quantile(EnsembleKalmanOutput* current_state,
                                                                    const Parameters &conditioned_on_parameters)
{
  return arma::colvec();
}
*/

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
  /*
  std::vector<double> new_schedule;
  new_schedule.reserve(2);
  new_schedule.push_back(this->current_value);
  new_schedule.push_back(parameters_in[this->variable_name][0]);
  this->schedule = new_schedule;
  this->mileometer.current_index.back() = 0;
  this->mileometer.limits.back() = 2;
  
  this->mileometer.increment();
  */
  
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
