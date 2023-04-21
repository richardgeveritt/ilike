#include "always_smc_termination.h"

AlwaysSMCTermination::AlwaysSMCTermination()
  :SMCTermination()
{
}

AlwaysSMCTermination::~AlwaysSMCTermination()
{
  
}

//Copy constructor for the AlwaysSMCTermination class.
AlwaysSMCTermination::AlwaysSMCTermination(const AlwaysSMCTermination &another)
  :SMCTermination(another)
{
  this->make_copy(another);
}

void AlwaysSMCTermination::operator=(const AlwaysSMCTermination &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  SMCTermination::operator=(another);
  this->make_copy(another);
}

SMCTermination* AlwaysSMCTermination::duplicate(void)const
{
  return( new AlwaysSMCTermination(*this));
}

void AlwaysSMCTermination::make_copy(const AlwaysSMCTermination &another)
{
}

bool AlwaysSMCTermination::terminate(double score)
{
  return true;
}
