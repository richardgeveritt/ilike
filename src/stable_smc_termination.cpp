#include "stable_smc_termination.h"
#include "utils.h"

StableSMCTermination::StableSMCTermination()
  :SMCTermination()
{
}

StableSMCTermination::StableSMCTermination(size_t number_in_a_row_in,
                                           double threshold_in)
:SMCTermination()
{
  this->number_in_a_row = number_in_a_row_in;
  this->threshold = threshold_in;
  this->counter = 0;
}

StableSMCTermination::~StableSMCTermination()
{
  
}

//Copy constructor for the StableSMCTermination class.
StableSMCTermination::StableSMCTermination(const StableSMCTermination &another)
  :SMCTermination(another)
{
  this->make_copy(another);
}

void StableSMCTermination::operator=(const StableSMCTermination &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCTermination::operator=(another);
  this->make_copy(another);
}

SMCTermination* StableSMCTermination::duplicate() const
{
  return( new StableSMCTermination(*this));
}

void StableSMCTermination::make_copy(const StableSMCTermination &another)
{
  this->counter = another.counter;
  this->threshold = another.threshold;
  this->number_in_a_row = another.number_in_a_row;
}

bool StableSMCTermination::terminate(double score)
{
  if (abs(score)<=this->threshold)
    this->counter = this->counter + 1;
  else
    this->counter = 0;
  
  if (this->counter >= this->number_in_a_row)
    return TRUE;
  else
    return FALSE;
}
