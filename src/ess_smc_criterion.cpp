#include "ess_smc_criterion.h"
#include "utils.h"
#include "smc_output.h"
#include "ensemble_kalman_output.h"

ESSSMCCriterion::ESSSMCCriterion()
  :SMCCriterion()
{
}

ESSSMCCriterion::ESSSMCCriterion(double desired_criterion_in)
  :SMCCriterion(desired_criterion_in)
{
  this->desired_criterion = desired_criterion_in;
}

ESSSMCCriterion::~ESSSMCCriterion()
{
  
}

//Copy constructor for the ESSSMCCriterion class.
ESSSMCCriterion::ESSSMCCriterion(const ESSSMCCriterion &another)
  :SMCCriterion(another)
{
  this->make_copy(another);
}

void ESSSMCCriterion::operator=(const ESSSMCCriterion &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCCriterion::operator=(another);
  this->make_copy(another);
}

SMCCriterion* ESSSMCCriterion::duplicate() const
{
  return( new ESSSMCCriterion(*this));
}

void ESSSMCCriterion::make_copy(const ESSSMCCriterion &another)
{
}

double ESSSMCCriterion::operator()(const Particles &particles) const
{
  if ( sum(particles.unnormalised_log_weights==-arma::datum::inf)==particles.unnormalised_log_weights.size() )
  {
    return(-this->desired_criterion);
  }
  else if ( sum(particles.unnormalised_log_weights==arma::datum::inf)==particles.unnormalised_log_weights.size() )
  {
    return(double(particles.incremental_log_weights.size())-this->desired_criterion);
  }
  else
  {
    return(particles.ess - this->desired_criterion);
  }
}

void ESSSMCCriterion::subsample_find_desired_criterion(SMCOutput* current_state)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

/*
void ESSSMCCriterion::find_desired_criterion(SMCOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

void ESSSMCCriterion::subsample_find_desired_criterion(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}
*/

void ESSSMCCriterion::find_desired_criterion(SMCOutput* current_state)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

double ESSSMCCriterion::operator()(const Ensemble &particles) const
{
  if ( sum(particles.unnormalised_log_weights==-arma::datum::inf)==particles.unnormalised_log_weights.size() )
  {
    return(-this->desired_criterion);
  }
  else if ( sum(particles.unnormalised_log_weights==arma::datum::inf)==particles.unnormalised_log_weights.size() )
  {
    return(double(particles.unnormalised_log_weights.size())-this->desired_criterion);
  }
  else
  {
    double ess = exp(2.0*log_sum_exp(particles.unnormalised_log_weights) - log_sum_exp(2.0*particles.unnormalised_log_weights));
    return(ess - this->desired_criterion);
  }
}

/*
void ESSSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}
*/

void ESSSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

bool ESSSMCCriterion::always_positive() const
{
  return false;
}
