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
          SMCCriterion* criterion_in,
          SMCTermination* termination_in)
{
  this->setup(the_worker_in,
              schedules_in,
              variable_names_in,
              criterion_in,
              termination_in);
}

void Sequencer::setup(SMCWorker* the_worker_in,
                      const std::vector< std::vector<double> > &schedules_in,
                      const std::vector<std::string> &variable_names_in,
                      SMCCriterion* criterion_in,
                      SMCTermination* termination_in)
{
  this->the_worker = the_worker_in;
  this->schedules = schedules_in;
  this->variable_names = variable_names_in;
  this->criterion = criterion_in;
  this->termination = termination_in;
  
  if ( (this->schedules.size()<1) || (this->schedules.size()<1) || (this->schedules.size()!=this->variable_names.size()) )
      Rcpp::stop("Sequencer::setup - invalid schedule.");
      
  if (this->schedules.back().back()>=this->schedules.back().front())
  {
    this->direction = 1.0;
  }
  else
  {
    this->direction = -1.0;
  }
  
  this->extra_bit = 0.001*(this->schedules.back().back()-this->schedules.back().front());
  
  std::vector<size_t> sizes;
  sizes.reserve(this->schedules.size());
  this->current_values.clear();
  this->current_values.reserve(this->schedules.size());
  for (std::vector< std::vector<double> >::const_iterator i=this->schedules.begin();
       i!=this->schedules.end();
       ++i)
  {
    if (i->size()<2)
      Rcpp::stop("Sequencer::setup - invalid schedule.");
    
    this->current_values.push_back(i->front());
    sizes.push_back(i->size());
  }
  
  this->mileometer = Mileometer(sizes);
  this->mileometer.increment();
  
  // Use epsilon_doubling in abc.r if we need help finding a maximum value to begin with.
}

Sequencer::Sequencer(const Sequencer &another)
{
  this->make_copy(another);
}

void Sequencer::operator=(const Sequencer &another)
{
  if(this == &another)
    return;
  
  this->schedules.clear();
  this->variable_names.clear();
  this->current_values.clear();
  
  if (this->criterion!=NULL)
    delete this->criterion;
  if (this->termination!=NULL)
    delete this->termination;

  this->make_copy(another);
}

void Sequencer::make_copy(const Sequencer &another)
{
  this->schedules = another.schedules;
  this->variable_names = another.variable_names;
  //this->criterion = another.criterion;
  //this->termination = another.termination;
  this->direction = another.direction;
  this->current_values = another.current_values;
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
}

void Sequencer::find_desired_criterion(SMCOutput* current_state)
{
  this->criterion->find_desired_criterion(current_state);
}

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

arma::colvec Sequencer::find_next_target_bisection(SMCOutput* current_state,
                                                   const Index* index)
{
  Parameters all_parameters = Parameters();
  for (size_t i=0; i<this->current_values.size(); ++i)
  {
    all_parameters[this->variable_names[i]] = this->current_values[i];
  }
  
  double current_bisect_value = this->current_values.back();
  std::vector<double> target_values = this->mileometer.get_current_values(this->schedules);
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    all_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->weight_for_adapting_sequence(current_state->back(),all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
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
    all_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->smcadaptive_given_smcfixed_weight(index,
                                                        current_state->back(),
                                                        all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->mileometer.increment();
  }
  
  return incremental_log_weights;
}

arma::colvec Sequencer::find_next_target_quantile(SMCOutput* current_state)
{
  return arma::colvec();
}

void Sequencer::find_desired_criterion(SMCOutput* current_state,
                                       const Parameters &conditioned_on_parameters)
{
  this->criterion->find_desired_criterion(current_state, conditioned_on_parameters);
}

// new method:
// sets desried ess via looking at current and last particle sets (based on \pi_t/\pi_t+1)
// does bisection with this target
// should also work for generic?

arma::colvec Sequencer::find_next_target_bisection(SMCOutput* current_state,
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
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;

  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;

    bisect_size = bisect_size/2.0;
    
    all_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->weight_for_adapting_sequence(current_state->back(),all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();

    this->current_score = (*this->criterion)(current_state->back());

    current_direction = std::copysign(1.0,this->current_score)*this->direction;

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
    all_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->smcadaptive_given_smcfixed_weight(index,
                                                        current_state->back(),
                                                        all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->mileometer.increment();
  }
  
  return incremental_log_weights;
}

arma::colvec Sequencer::find_next_target_quantile(SMCOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  return arma::colvec();
}

arma::colvec Sequencer::subsample_find_next_target_bisection(SMCOutput* current_state,
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
  double next_value = target_values.back()+this->extra_bit;
  double bisect_size = abs(next_value-current_bisect_value)/2.0;
  double current_direction = this->direction;
  double new_bisect_value;
  arma::colvec incremental_log_weights;
  
  for (size_t i=0; i<100; ++i)
  {
    new_bisect_value = current_bisect_value + current_direction*bisect_size;
    
    bisect_size = bisect_size/2.0;
    
    all_parameters[this->variable_names.back()] = new_bisect_value;
    
    // Call SMC specific weight here (done for MCMC, just target for PMC).
    this->the_worker->the_smc->subsample_weight_for_adapting_sequence(current_state->back(),all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->current_score = (*this->criterion)(current_state->back());
    
    current_direction = std::copysign(1.0,this->current_score)*this->direction;
    
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
    all_parameters[this->variable_names.back()] = new_bisect_value;
    this->the_worker->subsample_smcadaptive_given_smcfixed_weight(index,
                                                                  current_state->back(),
                                                                  all_parameters);
    
    incremental_log_weights = this->the_worker->get_unnormalised_log_incremental_weights();
    
    this->mileometer.increment();
  }
  
  return incremental_log_weights;
}

arma::colvec Sequencer::subsample_find_next_target_quantile(SMCOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  return arma::colvec();
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
  new_schedule.push_back(this->current_values.back());
  new_schedule.push_back(parameters_in[this->variable_names.back()][0]);
  this->schedules.back() = new_schedule;
  this->mileometer.current_index.back() = 0;
  this->mileometer.limits.back() = 2;
  
  this->mileometer.increment();
}
