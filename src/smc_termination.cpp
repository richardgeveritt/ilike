#include "smc_termination.h"

SMCTermination::SMCTermination()
{
}

SMCTermination::~SMCTermination()
{

}

SMCTermination::SMCTermination(const SMCTermination &another)
{
  this->make_copy(another);
}

void SMCTermination::operator=(const SMCTermination &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMCTermination::make_copy(const SMCTermination &another)
{
}
