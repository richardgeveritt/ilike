#include "sequencer.h"

Sequencer::Sequencer(void)
{
}

Sequencer::Sequencer(const Sequencer &another)
{
  this->make_copy(another);
}

Sequencer::~Sequencer(void)
{
}

void Sequencer::operator=(const Sequencer &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Sequencer::make_copy(const Sequencer &another)
{
  //Does nothing since no member variables to copy.
}
