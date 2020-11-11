#include "abc_sequencer.h"

//Default constructor.
ABCSequencer::ABCSequencer(void)
  :Sequencer()
{
}

//Copy constructor for the ABCSequencer class.
ABCSequencer::ABCSequencer(const ABCSequencer &another)
  :Sequencer(another)
{
  this->make_copy(another);
}

//Destructor for the ABCSequencer class.
ABCSequencer::~ABCSequencer(void)
{
}

void ABCSequencer::operator=(const ABCSequencer &another)
{
  if(this == &another){ //if a==a
    return;
  }

  Sequencer::operator=(another);
  this->make_copy(another);
}

Sequencer* ABCSequencer::duplicate(void)const
{
  return( new ABCSequencer(*this));
}

void ABCSequencer::make_copy(const ABCSequencer &another)
{

}
