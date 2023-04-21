#include "se_mcmc_termination.h"

SEMCMCTermination::SEMCMCTermination()
  :MCMCTermination()
{
}

SEMCMCTermination::SEMCMCTermination(double threshold_in,
                                     size_t max_number_of_iterations_in)
  :MCMCTermination()
{
  this->threshold = threshold_in;
  this->max_number_of_iterations = max_number_of_iterations_in;
}

SEMCMCTermination::~SEMCMCTermination()
{
  
}

//Copy constructor for the SEMCMCTermination class.
SEMCMCTermination::SEMCMCTermination(const SEMCMCTermination &another)
  :MCMCTermination(another)
{
  this->make_copy(another);
}

void SEMCMCTermination::operator=(const SEMCMCTermination &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MCMCTermination::operator=(another);
  this->make_copy(another);
}

MCMCTermination* SEMCMCTermination::duplicate(void)const
{
  return( new SEMCMCTermination(*this));
}

void SEMCMCTermination::make_copy(const SEMCMCTermination &another)
{
  this->threshold = another.threshold;
  this->max_number_of_iterations = another.max_number_of_iterations;
}

bool SEMCMCTermination::terminate()
{
  return true;
}
