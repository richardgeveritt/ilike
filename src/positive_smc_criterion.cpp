#include "positive_smc_criterion.h"

PositiveSMCCriterion::PositiveSMCCriterion()
  :SMCCriterion()
{
}

PositiveSMCCriterion::~PositiveSMCCriterion()
{
  
}

// Copy constructor for the PositiveSMCCriterion class.
PositiveSMCCriterion::PositiveSMCCriterion(const PositiveSMCCriterion &another)
  :SMCCriterion(another)
{
  this->make_copy(another);
}

void PositiveSMCCriterion::operator=(const PositiveSMCCriterion &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCCriterion::operator=(another);
  this->make_copy(another);
}

SMCCriterion* PositiveSMCCriterion::duplicate() const
{
  return( new PositiveSMCCriterion(*this));
}

void PositiveSMCCriterion::make_copy(const PositiveSMCCriterion &another)
{
}

double PositiveSMCCriterion::operator()(const Particles &particles) const
{
  return 1.0;
}

void PositiveSMCCriterion::find_desired_criterion(SMCOutput* current_state)
{
}

/*
void PositiveSMCCriterion::find_desired_criterion(SMCOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
}
*/

void PositiveSMCCriterion::subsample_find_desired_criterion(SMCOutput* current_state)
{
}

double PositiveSMCCriterion::operator()(const Ensemble &particles) const
{
  return 1.0;
}

void PositiveSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state)
{
  
}

/*
void PositiveSMCCriterion::find_desired_criterion(EnsembleKalmanOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  
}
*/

bool PositiveSMCCriterion::always_positive() const
{
  return true;
}
