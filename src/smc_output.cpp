#include "smc_output.h"

SMCOutput::SMCOutput()
{
}

SMCOutput::~SMCOutput()
{

}

//Copy constructor for the SMCOutput class.
SMCOutput::SMCOutput(const SMCOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void SMCOutput::operator=(const SMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* SMCOutput::duplicate(void)const
{
  return( new SMCOutput(*this));
}

void SMCOutput::make_copy(const SMCOutput &another)
{
  this->all_particles = another.all_particles;
  this->all_proposed = another.all_proposed;
  this->log_normalised_weights = another.log_normalised_weights;
  this->log_normalising_constant = another.log_normalising_constant;
}
