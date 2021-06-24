#include "smc_mcmc_move.h"

///Default constructor, Shape and Scale parameters have no values.
SMCMCMCMove::SMCMCMCMove(const ModelAndAlgorithm* model_and_algorithm_in)
  :SMC(model_and_algorithm_in)
{
}

//Copy constructor for the SMCMCMCMove class.
SMCMCMCMove::SMCMCMCMove(const SMCMCMCMove &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCMCMCMove class.
SMCMCMCMove::~SMCMCMCMove(void)
{
}

void SMCMCMCMove::operator=(const SMCMCMCMove &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCMCMCMove::duplicate(void) const
{
  return( new SMCMCMCMove(*this));
}

void SMCMCMCMove::make_copy(const SMCMCMCMove &another)
{

}

void SMCMCMCMove::smc_step(void)
{

}

void SMCMCMCMove::weight_update(void)
{

}
