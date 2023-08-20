#include "sequencer.h"
#include "smc_criterion.h"
#include "smc_termination.h"
#include "smc_worker.h"
#include "smc_output.h"
#include "smc.h"
#include "index.h"

Sequencer::Sequencer()
{
  this->criterion = NULL;
  this->termination = NULL;
  this->number_of_bisections = 25;
}

Sequencer::~Sequencer()
{
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;
}

//Sequencer::Sequencer(const std::vector<double> &schedule_in,
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

Sequencer::Sequencer(SMCWorker* the_worker_in,
                     const std::vector<double> &schedule_in,
                     const std::string &variable_in,
                     size_t number_of_bisections_in,
                     SMCCriterion* criterion_in,
                     SMCTermination* termination_in)
{
  std::vector< std::vector<double> > schedules_in;
  schedules_in.push_back(schedule_in);
  std::vector<std::string> variable_names_in;
  variable_names_in.push_back(variable_in);
  this->setup(the_worker_in,
              schedules_in,
              variable_names_in,
              number_of_bisections_in,
              criterion_in,
              termination_in);
}

//Sequencer::Sequencer(const std::vector< std::vector<double> > &schedules_in,
//          const std::vector<std::string> &variable_names_in)
//{
//  SMCCriterion* criterion_in = NULL;
//  double desired_criterion_in = 0.0;
//  this->setup(schedules_in, variable_names_in, desired_criterion_in, criterion_in);
//}

Sequencer::Sequencer(SMCWorker* the_worker_in,
                     const std::vector< std::vector<double> > &schedules_in,
                     const std::vector<std::string> &variable_names_in,
                     size_t number_of_bisections_in,
                     SMCCriterion* criterion_in,
                     SMCTermination* termination_in)
{
  this->setup(the_worker_in,
              schedules_in,
              variable_names_in,
              number_of_bisections_in,
              criterion_in,
              termination_in);
}

void Sequencer::setup(SMCWorker* the_worker_in,
                      const std::vector< std::vector<double> > &schedules_in,
                      const std::vector<std::string> &variable_names_in,
                      size_t number_of_bisections_in,
                      SMCCriterion* criterion_in,
                      SMCTermination* termination_in)
{
  this->the_worker = the_worker_in;
  this->schedules = schedules_in;
  this->variable_names = variable_names_in;
  this->criterion = criterion_in;
  this->termination = termination_in;
  
  this->number_of_bisections = number_of_bisections_in;
  
  if ( (this->schedules.size()<1) || (this->schedules.size()<1) || (this->schedules.size()!=this->variable_names.size()) )
      Rcpp::stop("Sequencer::setup - invalid schedule.");
  
  std::vector<size_t> sizes;
  sizes.reserve(this->schedules.size());
  use_final.reserve(this->schedules.size());
  for (size_t i=0;
       i<this->schedules.size();
       ++i)
  {
    if (this->schedules[i].size()<1)
      Rcpp::stop("Sequencer::setup - invalid schedule.");
    sizes.push_back(this->schedules[i].size());
    
    use_final.push_back(true);
  }
  this->mileometer = Mileometer(sizes);
      
  this->reset();
}

void Sequencer::set_initial_schedule_parameters()
{
  this->mileometer.reset();
  this->mileometer.increment();
  
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    double value_to_use;
    if ((this->use_final[i]==true) && (this->mileometer[i]==0))
    {
      value_to_use = this->schedules[i][this->schedules[i].size()-1];
    }
    else
    {
      value_to_use = this->schedules[i][this->mileometer[i]];
    }
    this->schedule_parameters[this->variable_names[i]] = value_to_use;
    
    if (i==this->variable_names.size()-1)
    {
      // do I need use_final?
      this->current_bisect_value = this->schedules[i][this->mileometer[i]];
    }
  }
}

void Sequencer::set_schedule_parameters()
{
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    double value_to_use;
    if ((this->use_final[i]==true) && (this->mileometer[i]==0))
    {
      value_to_use = this->schedules[i][this->schedules[i].size()-1];
    }
    else
    {
      value_to_use = this->schedules[i][this->mileometer[i]];
    }
    this->schedule_parameters[this->variable_names[i]] = value_to_use;
    
    this->mileometer.increment();
    
    if (i==this->variable_names.size()-1)
    {
      // do I need use_final?
      this->current_bisect_value = this->schedules[i][this->mileometer[i]];
    }
  }
}

Sequencer::Sequencer(const Sequencer &another)
{
  this->make_copy(another);
}

Sequencer& Sequencer::operator=(const Sequencer &another)
{
  if(this == &another)
    return *this;
  
  this->schedules.clear();
  this->variable_names.clear();
  //this->current_values.clear();
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;

  this->make_copy(another);
  
  return *this;
}

void Sequencer::make_copy(const Sequencer &another)
{
  this->schedules = another.schedules;
  this->variable_names = another.variable_names;
  this->use_final = another.use_final;
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = another.direction;
  //this->current_values = another.current_values;
  this->mileometer = another.mileometer;
  this->current_score = another.current_score;
  this->current_bisect_value = another.current_bisect_value;
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
  this->number_of_bisections = another.number_of_bisections;
  this->schedule_parameters = another.schedule_parameters;
}

Sequencer::Sequencer(Sequencer &&another)
{
  this->make_copy(std::move(another));
}

Sequencer& Sequencer::operator=(Sequencer &&another)
{
  if(this == &another)
    return *this;
  
  this->schedules.clear();
  this->variable_names.clear();
  //this->current_values.clear();
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;
  
  this->make_copy(std::move(another));
  
  return *this;
}

void Sequencer::make_copy(Sequencer &&another)
{
  this->schedules = std::move(another.schedules);
  this->variable_names = std::move(another.variable_names);
  this->use_final = std::move(another.use_final);
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = std::move(another.direction);
  //this->current_values = std::move(another.current_values);
  this->mileometer = std::move(another.mileometer);
  this->current_score = std::move(another.current_score);
  this->current_bisect_value = std::move(another.current_bisect_value);
  this->number_of_bisections = another.number_of_bisections;
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
  this->schedule_parameters = std::move(another.schedule_parameters);
                               
  another.schedules = std::vector< std::vector<double> >();
  another.variable_names = std::vector<std::string>();
  another.direction = 0.0;
  //another.current_values = std::vector<double>();
  another.mileometer = Mileometer();
  another.current_bisect_value = 0.0;
  another.current_score = 0.0;
  another.the_worker = NULL;
  another.criterion = NULL;
  another.termination = NULL;
  another.schedule_parameters = another.schedule_parameters;
  another.number_of_bisections = 0;
}

void Sequencer::find_desired_criterion(SMCOutput* current_state)
{
  this->criterion->find_desired_criterion(current_state);
}

void Sequencer::subsample_find_desired_criterion(SMCOutput* current_state)
{
  this->criterion->subsample_find_desired_criterion(current_state);
}

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

void Sequencer::find_next_target_bisection(SMCOutput* current_state,
                                           const Index* index)
{
  //Parameters all_parameters = Parameters();
  //for (size_t i=0; i<this->current_values.size(); ++i)
  //{
  //  all_parameters[this->variable_names[i]] = this->current_values[i];
  //}
  
  //this->set_schedule_parameters();
  
  //double current_bisect_value = this->current_values.back();
  //std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  //arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  //this->schedule_parameters[this->variable_names.back()] = target_values.back();
  
  // if we want to make this generic, need to hide the next two lines in a function that can have different choices
  
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->current_bisect_value>=this->schedules.back()[this->mileometer.back()])
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    if (this->current_score>=0.0)
    {
      //this->mileometer.increment();
      this->set_schedule_parameters();
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
    
    // if there is only one value in the schedule, but we did not get to the end in one go, then we need to determine the start point in order to do the bisection.
    
    // first identify if there is only one value in the schedule
    if (this->schedules.back().size()==1)
    {
      // we will assume that we need to find an upper bound
      // perform a doubling routine, which we keep doing until we obtain an upper bound (max 300 iterations)
      
      //std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
      //current_bisect_value = target_values[0];
      if (this->current_bisect_value==0.0)
        this->current_bisect_value = this->current_bisect_value + 0.00001;
      double new_bisect_value = 0.0;
      
      for (size_t i=0; i<300; ++i)
      {
        new_bisect_value = 2.0*this->current_bisect_value;
        
        this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
        
        // Call SMC specific weight here (done for MCMC, just target for PMC).
        this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                                current_state->back());
        current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
        
        this->current_score = (*this->criterion)(current_state->back());
        
        if (this->current_score>=0.0)
        {
          break;
        }
        
        this->current_bisect_value = new_bisect_value;
      }
      
      // change sequencer
      std::vector<double> new_final_values;
      new_final_values.push_back(new_bisect_value);
      new_final_values.push_back(this->schedules.back()[0]);
      this->schedules.back() = new_final_values;
      
      this->mileometer.reset_final_dimension(2);
      if (new_final_values[1]>new_final_values[0])
      {
        this->direction = 1.0;
      }
      else
      {
        this->direction = -1.0;
      }
      //this->current_values.back() = new_bisect_value;
      
      //this->mileometer.increment();
      this->set_schedule_parameters();
      
      //current_bisect_value = this->current_values.back();
      //target_values = this->mileometer.get_current_values(this->schedules);
    }
    
    //next_value = this->current_bisect_value;//schedule_parameters[this->variable_names.back()][0];
    this->current_bisect_value = this->schedules.back()[this->mileometer.back()-1];
  }
  else
  {
    this->schedule_parameters[this->variable_names.back()] = this->schedules.back()[this->mileometer.back()];
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    if (this->current_score>=0.0)
    {
      //this->mileometer.increment();
      this->set_schedule_parameters();
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
  }
  
  double next_value = this->schedules.back()[this->mileometer.back()];
  // target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  //this->current_values.back() = new_bisect_value;
  
  // Check if we went past the end.
  /*
  if ((std::copysign(1.0,this->current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
  }
  */
  
  //return incremental_log_weights;
}

void Sequencer::find_next_target_quantile(SMCOutput* current_state,
                                          const Index* index)
{
  // for now, just move on to the next target in the list
  /*
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  this->schedule_parameters[this->variable_names.back()] = target_values.back();
  
  // if we want to make this generic, need to hide the next two lines in a function that can have different choices
  this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                          current_state->back());
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->back());
  if (this->current_score>=0.0)
  {
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = target_values.back();
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    return;
  }
  else
  {
    Rcpp::stop("Sequencer::find_next_target_quantile - not yet implemented.");
  }
  */
  
  Rcpp::stop("Sequencer::find_next_target_quantile - not yet implemented.");
}

/*
void Sequencer::find_desired_criterion(SMCOutput* current_state,
                                       const Parameters &conditioned_on_parameters)
{
  this->criterion->find_desired_criterion(current_state, conditioned_on_parameters);
}

void Sequencer::subsample_find_desired_criterion(SMCOutput* current_state,
                                                 const Parameters &conditioned_on_parameters)
{
  this->criterion->subsample_find_desired_criterion(current_state, conditioned_on_parameters);
}
*/

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

/*
void Sequencer::find_next_target_bisection(SMCOutput* current_state,
                                           const Index* index,
                                           const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = conditioned_on_parameters;
  for (size_t i=0; i<this->current_values.size(); ++i)
  {
    all_parameters[this->variable_names[i]] = this->current_values[i];
  }
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  all_parameters[this->variable_names.back()] = target_values.back();
  this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                          current_state->back(),
                                                          all_parameters);
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->back());
  if (this->current_score>0.0)
  {
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = target_values.back();
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    return;
  }
  
  // if there is only one value in the schedule, but we did not get to the end in one go, then we need to determine the start point in order to do the bisection.
  
  // first identify if there is only one value in the schedule
  if (this->schedules.back().size()==1)
  {
    // we will assume that we need to find an upper bound
    // perform a doubling routine, which we keep doing until we obtain an upper bound (max 300 iterations)
    
    std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
    current_bisect_value = target_values[0];
    if (current_bisect_value==0.0)
      current_bisect_value = current_bisect_value + 0.00001;
    double new_bisect_value = 0.0;
    
    for (size_t i=0; i<300; ++i)
    {
      new_bisect_value = 2.0*current_bisect_value;
      
      all_parameters[this->variable_names.back()] = new_bisect_value;
      
      // Call SMC specific weight here (done for MCMC, just target for PMC).
      this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                              current_state->back(),
                                                              all_parameters);
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
      
      this->current_score = (*this->criterion)(current_state->back());
      
      if (this->current_score>=0.0)
      {
        break;
      }
      
      current_bisect_value = new_bisect_value;
    }
    
    // change sequencer
    std::vector<double> new_final_values;
    new_final_values.push_back(new_bisect_value);
    new_final_values.push_back(target_values[0]);
    this->schedules.back() = new_final_values;
    
    this->mileometer.reset_final_dimension(2);
    if (new_final_values[1]>new_final_values[0])
    {
      this->direction = 1.0;
    }
    else
    {
      this->direction = -1.0;
    }
    this->current_values.back() = new_bisect_value;
    
    this->mileometer.increment();
    
    current_bisect_value = this->current_values.back();
    target_values = this->mileometer.get_current_values(this->schedules);
  }
  
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    all_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back(),
                                                            all_parameters);
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_values.back() = new_bisect_value;
  this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    all_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                            current_state->back(),
                                                            all_parameters);
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
  }
  
  //return incremental_log_weights;
}
*/

/*
void Sequencer::find_next_target_quantile(SMCOutput* current_state,
                                          const Index* index)
{
  Parameters all_parameters = conditioned_on_parameters;
  for (size_t i=0; i<this->current_values.size(); ++i)
  {
    all_parameters[this->variable_names[i]] = this->current_values[i];
  }
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  all_parameters[this->variable_names.back()] = target_values.back();
  this->the_worker->the_smc->weight_for_adapting_sequence(index,
                                                          current_state->back(),
                                                          all_parameters);
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->back());
  if (this->current_score>0.0)
  {
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = target_values.back();
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    return;
  }
  else
  {
    Rcpp::stop("Sequencer::find_next_target_quantile - not yet implemented.");
  }
}
*/

void Sequencer::subsample_find_next_target_bisection(SMCOutput* current_state,
                                                     const Index* index)
{
  //arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  //this->schedule_parameters[this->variable_names.back()] = target_values.back();
  
  // if we want to make this generic, need to hide the next two lines in a function that can have different choices
  
  // if we have not just found a value for the parameter that is at one of the points in the schedule
  if (this->current_bisect_value>=this->schedules.back()[this->mileometer.back()])
  {
    // check to see if we already reached the next point in the schedule
    
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                      current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    if (this->current_score>=0.0)
    {
      this->set_schedule_parameters();
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
    
    // if there is only one value in the schedule, but we did not get to the end in one go, then we need to determine the start point in order to do the bisection.
    
    // first identify if there is only one value in the schedule
    if (this->schedules.back().size()==1)
    {
      // we will assume that we need to find an upper bound
      // perform a doubling routine, which we keep doing until we obtain an upper bound (max 300 iterations)
      
      //std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
      //current_bisect_value = target_values[0];
      if (this->current_bisect_value==0.0)
        this->current_bisect_value = this->current_bisect_value + 0.00001;
      double new_bisect_value = 0.0;
      
      for (size_t i=0; i<300; ++i)
      {
        new_bisect_value = 2.0*this->current_bisect_value;
        
        this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
        
        // Call SMC specific weight here (done for MCMC, just target for PMC).
        this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                          current_state->back());
        current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
        
        this->current_score = (*this->criterion)(current_state->back());
        
        if (this->current_score>=0.0)
        {
          break;
        }
        
        this->current_bisect_value = new_bisect_value;
      }
      
      // change sequencer
      std::vector<double> new_final_values;
      new_final_values.push_back(new_bisect_value);
      new_final_values.push_back(this->schedules.back()[0]);
      this->schedules.back() = new_final_values;
      
      this->mileometer.reset_final_dimension(2);
      if (new_final_values[1]>new_final_values[0])
      {
        this->direction = 1.0;
      }
      else
      {
        this->direction = -1.0;
      }
      //this->current_values.back() = new_bisect_value;
      
      this->mileometer.increment();
      this->set_schedule_parameters();
      
      //current_bisect_value = this->current_values.back();
      //target_values = this->mileometer.get_current_values(this->schedules);
    }
    
    //next_value = this->current_bisect_value;//schedule_parameters[this->variable_names.back()][0];
    this->current_bisect_value = this->schedules.back()[this->mileometer.back()-1];
  }
  else
  {
    this->schedule_parameters[this->variable_names.back()] = this->schedules.back()[this->mileometer.back()];
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                      current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    if (this->current_score>=0.0)
    {
      this->mileometer.increment();
      this->set_schedule_parameters();
      //this->schedule_parameters[this->variable_names.back()] = target_values.back();
      
      //this->mileometer.increment();
      //this->current_values = this->mileometer.get_current_values(this->schedules);
      
      return;
    }
  }
  
  double next_value = this->schedules.back()[this->mileometer.back()];
  // target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-this->current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<this->number_of_bisections; ++i)
  {
    new_bisect_value = this->current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                      current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    this->current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  /*
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  this->schedule_parameters[this->variable_names.back()] = target_values.back();
  
  
  this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                    current_state->back());
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->back());
  if (this->current_score>0.0)
  {
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = target_values.back();
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    return;
  }
  
  // if there is only one value in the schedule, but we did not get to the end in one go, then we need to determine the start point in order to do the bisection.
  
  // first identify if there is only one value in the schedule
  if (this->schedules.back().size()==1)
  {
    // we will assume that we need to find an upper bound
    // perform a doubling routine, which we keep doing until we obtain an upper bound (max 300 iterations)
    
    std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
    current_bisect_value = target_values[0];
    if (current_bisect_value==0.0)
      current_bisect_value = current_bisect_value + 0.00001;
    double new_bisect_value = 0.0;
    
    for (size_t i=0; i<300; ++i)
    {
      new_bisect_value = 2.0*current_bisect_value;
      
      this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
      
      // Call SMC specific weight here (done for MCMC, just target for PMC).
      this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                        current_state->back());
      current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
      
      this->current_score = (*this->criterion)(current_state->back());
      
      if (this->current_score>=0.0)
      {
        break;
      }
      
      current_bisect_value = new_bisect_value;
    }
    
    // change sequencer
    std::vector<double> new_final_values;
    new_final_values.push_back(new_bisect_value);
    new_final_values.push_back(target_values[0]);
    this->schedules.back() = new_final_values;
    
    this->mileometer.reset_final_dimension(2);
    if (new_final_values[1]>new_final_values[0])
    {
      this->direction = 1.0;
    }
    else
    {
      this->direction = -1.0;
    }
    this->current_values.back() = new_bisect_value;
    
    this->mileometer.increment();
    
    current_bisect_value = this->current_values.back();
    target_values = this->mileometer.get_current_values(this->schedules);
  }
  
  double next_value = target_values.back();//+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value = 0.0;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                      current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
    current_bisect_value = new_bisect_value;
    
    if (this->current_score==0.0)
    {
      break;
    }
    
  }
  
  this->current_values.back() = new_bisect_value;
  this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
  
  // Check if we went past the end.
  if ((std::copysign(1.0,current_bisect_value-target_values.back())*this->direction)>0)
  {
    new_bisect_value = target_values.back();
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                      current_state->back());
    current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = new_bisect_value;
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
  }
  
  //return incremental_log_weights;
  */
}

void Sequencer::subsample_find_next_target_quantile(SMCOutput* current_state,
                                                    const Index* index)
{
  /*
  this->set_schedule_parameters();
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  arma::colvec incremental_log_weights;
  
  // see if we need to do the bisection
  this->schedule_parameters[this->variable_names.back()] = target_values.back();
  
  this->the_worker->the_smc->subsample_weight_for_adapting_sequence(index,
                                                                    current_state->back());
  current_state->back().update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  this->current_score = (*this->criterion)(current_state->back());
  if (this->current_score>0.0)
  {
    //incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    this->set_schedule_parameters();
    this->schedule_parameters[this->variable_names.back()] = target_values.back();
    
    this->mileometer.increment();
    this->current_values = this->mileometer.get_current_values(this->schedules);
    
    return;
  }
  else
  {
    Rcpp::stop("Sequencer::find_next_target_quantile - not yet implemented.");
  }
  */
  
  Rcpp::stop("Sequencer::find_next_target_quantile - not yet implemented.");
}

bool Sequencer::check_termination()
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

void Sequencer::set_next_with_parameter(const Parameters &parameters_in)
{
  std::vector<double> new_schedule;
  new_schedule.reserve(2);
  new_schedule.push_back(this->schedules.back()[this->mileometer.back()]);
  new_schedule.push_back(parameters_in[this->variable_names.back()][0]);
  this->schedules.back() = new_schedule;
  this->mileometer.reset_final_dimension(2);
  if (this->schedules.back().back()>=this->schedules.back().front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }
  
  this->mileometer.increment();
}

void Sequencer::reset()
{
  if (this->schedules.back().back()>=this->schedules.back().front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }

  //this->mileometer.increment();
  this->set_initial_schedule_parameters();
  
  //this->extra_bit = 0.001*(this->schedules.back().back()-this->schedules.back().front());
  
  //std::vector<size_t> sizes;
  //sizes.reserve(this->schedules.size());
  //this->current_values.clear();
  //this->current_values.reserve(this->schedules.size());
  /*
  for (size_t i=0;
       i<this->schedules.size();
       ++i)
  {
    //if (i->size()<1)
    //  Rcpp::stop("Sequencer::setup - invalid schedule.");
    this->schedule_parameters[this->variable_names[i]] = this->schedules[i][0]
    
    this->current_values.push_back(i->front());
    //sizes.push_back(i->size());
  }
  */
  
  //this->mileometer = std::move(Mileometer(sizes));
  
  
  
  // Use epsilon_doubling in abc.r if we need help finding a maximum value to begin with.
}
