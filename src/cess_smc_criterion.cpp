#include "cess_smc_criterion.h"
#include "utils.h"
#include "smc_output.h"
#include "ensemble_kalman_output.h"

CESSSMCCriterion::CESSSMCCriterion()
  :SMCCriterion()
{
}

CESSSMCCriterion::CESSSMCCriterion(double desired_criterion_in)
  :SMCCriterion(desired_criterion_in)
{
  this->desired_criterion = desired_criterion_in;
}

CESSSMCCriterion::~CESSSMCCriterion()
{
  
}

//Copy constructor for the CESSSMCCriterion class.
CESSSMCCriterion::CESSSMCCriterion(const CESSSMCCriterion &another)
  :SMCCriterion(another)
{
  this->make_copy(another);
}

void CESSSMCCriterion::operator=(const CESSSMCCriterion &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCCriterion::operator=(another);
  this->make_copy(another);
}

SMCCriterion* CESSSMCCriterion::duplicate(void)const
{
  return( new CESSSMCCriterion(*this));
}

void CESSSMCCriterion::make_copy(const CESSSMCCriterion &another)
{
}

double CESSSMCCriterion::operator()(const Particles &particles) const
{
  if ( sum(particles.incremental_log_weights==R_NegInf)==particles.previous_normalised_log_weights.size() )
  {
    return(-this->desired_criterion);
  }
  else if ( sum(particles.incremental_log_weights==R_PosInf)==particles.previous_normalised_log_weights.size() )
  {
    return(double(particles.incremental_log_weights.size())-this->desired_criterion);
  }
  else
  {
    return(exp(log(particles.incremental_log_weights.size()) + 2.0*(log_sum_exp(particles.previous_normalised_log_weights + particles.incremental_log_weights)) - log_sum_exp(particles.previous_normalised_log_weights + 2.0*particles.incremental_log_weights)) - this->desired_criterion);
  }
}

void CESSSMCCriterion::find_desired_criterion(SMCOutput* current_state)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

void CESSSMCCriterion::find_desired_criterion(SMCOutput* current_state,
                                              const Parameters &conditioned_on_parameters)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

double CESSSMCCriterion::operator()(const Ensemble &particles) const
{
  throw std::runtime_error("ESSSMCCriterion::operator() - not written yet.");
}

void CESSSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}

void CESSSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state,
                                              const Parameters &conditioned_on_parameters)
{
  // use ratio method from https://arxiv.org/pdf/1907.01505.pdf
  
  // different at first step
  
  // use previous_particles.previous_llhd
  // use current_particles.current_llhd
}
