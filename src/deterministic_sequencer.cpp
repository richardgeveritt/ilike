#include "deterministic_sequencer.h"

///Default constructor, Shape and Scale parameters have no values.
DeterministicSequencer::DeterministicSequencer(void)
  :Sequencer()
{
}

//Copy constructor for the DeterministicSequencer class.
DeterministicSequencer::DeterministicSequencer(const DeterministicSequencer &another)
  :Sequencer(another)
{
  this->make_copy(another);
}

//Destructor for the DeterministicSequencer class.
DeterministicSequencer::~DeterministicSequencer(void)
{
}

void DeterministicSequencer::operator=(const DeterministicSequencer &another)
{
  if(this == &another){ //if a==a
    return;
  }

  Sequencer::operator=(another);
  this->make_copy(another);
}

Sequencer* DeterministicSequencer::duplicate(void)const
{
  return( new DeterministicSequencer(*this));
}

void DeterministicSequencer::make_copy(const DeterministicSequencer &another)
{

}
