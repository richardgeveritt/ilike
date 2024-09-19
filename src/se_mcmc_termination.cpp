#include "se_mcmc_termination.h"
#include "mcmc.h"

namespace ilike
{
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
  //this->counter = counter_pointer;
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

MCMCTermination* SEMCMCTermination::duplicate() const
{
  return( new SEMCMCTermination(*this));
}

void SEMCMCTermination::make_copy(const SEMCMCTermination &another)
{
  this->threshold = another.threshold;
  this->max_number_of_iterations = another.max_number_of_iterations;
  this->counter = another.counter;
}

bool SEMCMCTermination::terminate()
{
  if (*this->counter==this->max_number_of_iterations)
  {
    *this->counter = 0;
    return true;
  }
  else
  {
    return false;
  }
}

void SEMCMCTermination::set_parameters(StandardMCMCOutput* mcmc_output)
{
  this->counter = mcmc_output->get_iteration_counter_pointer();
}
}
