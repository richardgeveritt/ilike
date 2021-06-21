#include "sequential_batch_simulator.h"

//Default constructor.
SequentialBatchSimulator::SequentialBatchSimulator(void)
  :BatchSimulator()
{
}

//Copy constructor for the SequentialBatchSimulator class.
SequentialBatchSimulator::SequentialBatchSimulator(const SequentialBatchSimulator &another)
  :BatchSimulator(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialBatchSimulator class.
SequentialBatchSimulator::~SequentialBatchSimulator(void)
{
}

void SequentialBatchSimulator::operator=(const SequentialBatchSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }

  BatchSimulator::operator=(another);
  this->make_copy(another);
}

BatchSimulator* SequentialBatchSimulator::duplicate(void)const
{
  return( new SequentialBatchSimulator(*this));
}

void SequentialBatchSimulator::make_copy(const SequentialBatchSimulator &another)
{

}

Particles SequentialBatchSimulator::simulate(void) const
{
  return Particles();
}
