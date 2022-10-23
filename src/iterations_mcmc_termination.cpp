#include "iterations_mcmc_termination.h"

IterationsMCMCTermination::IterationsMCMCTermination()
  :MCMCTermination()
{
}

IterationsMCMCTermination::IterationsMCMCTermination(size_t number_of_iterations_in,
                                                     size_t* counter_pointer)
:MCMCTermination()
{
  this->number_of_iterations = number_of_iterations_in;
  this->counter = counter_pointer;
}

IterationsMCMCTermination::~IterationsMCMCTermination()
{
  
}

//Copy constructor for the IterationsMCMCTermination class.
IterationsMCMCTermination::IterationsMCMCTermination(const IterationsMCMCTermination &another)
  :MCMCTermination(another)
{
  this->make_copy(another);
}

void IterationsMCMCTermination::operator=(const IterationsMCMCTermination &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  MCMCTermination::operator=(another);
  this->make_copy(another);
}

MCMCTermination* IterationsMCMCTermination::duplicate(void)const
{
  return( new IterationsMCMCTermination(*this));
}

void IterationsMCMCTermination::make_copy(const IterationsMCMCTermination &another)
{
  this->number_of_iterations = another.number_of_iterations;
  this->counter = another.counter;
}

bool IterationsMCMCTermination::terminate()
{
  *this->counter = *this->counter + 1;
  if (*this->counter==this->number_of_iterations)
  {
    *this->counter = 0;
    return true;
  }
  else
  {
    return false;
  }
}
