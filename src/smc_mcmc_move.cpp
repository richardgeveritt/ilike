#include "smc_mcmc_move.h"

///Default constructor, Shape and Scale parameters have no values.
SMCMCMCMove::SMCMCMCMove(void)
  :SMC()
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

SMC* SMCMCMCMove::duplicate(void)const
{
  return( new SMCMCMCMove(*this));
}

void SMCMCMCMove::make_copy(const SMCMCMCMove &another)
{

}
