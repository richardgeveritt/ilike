#include "annealing_sequencer.h"

//Default constructor.
AnnealingSequencer::AnnealingSequencer(void)
  :Sequencer()
{
}

//Copy constructor for the AnnealingSequencer class.
AnnealingSequencer::AnnealingSequencer(const AnnealingSequencer &another)
  :Sequencer(another)
{
  this->make_copy(another);
}

//Destructor for the AnnealingSequencer class.
AnnealingSequencer::~AnnealingSequencer(void)
{
}

void AnnealingSequencer::operator=(const AnnealingSequencer &another)
{
  if(this == &another){ //if a==a
    return;
  }

  Sequencer::operator=(another);
  this->make_copy(another);
}

Sequencer* AnnealingSequencer::duplicate(void)const
{
  return( new AnnealingSequencer(*this));
}

void AnnealingSequencer::make_copy(const AnnealingSequencer &another)
{

}
