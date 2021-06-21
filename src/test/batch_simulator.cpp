#include "batch_simulator.h"

BatchSimulator::BatchSimulator(void)
{
}

BatchSimulator::BatchSimulator(const BatchSimulator &another)
{
  this->make_copy(another);
}

BatchSimulator::~BatchSimulator(void)
{
}

void BatchSimulator::operator=(const BatchSimulator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void BatchSimulator::make_copy(const BatchSimulator &another)
{
  //Does nothing since no member variables to copy.
}
