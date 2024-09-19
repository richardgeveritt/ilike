#include "mcmc_termination.h"

namespace ilike
{
MCMCTermination::MCMCTermination()
{
}

MCMCTermination::~MCMCTermination()
{
  
}

MCMCTermination::MCMCTermination(const MCMCTermination &another)
{
  this->make_copy(another);
}

void MCMCTermination::operator=(const MCMCTermination &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void MCMCTermination::make_copy(const MCMCTermination &another)
{
}
}
