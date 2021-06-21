#include "rcppparallel_batch_simulator.h"

//Default constructor.
RcppParallelBatchSimulator::RcppParallelBatchSimulator(void)
  :BatchSimulator()
{
}

//Copy constructor for the RcppParallelBatchSimulator class.
RcppParallelBatchSimulator::RcppParallelBatchSimulator(const RcppParallelBatchSimulator &another)
  :BatchSimulator(another)
{
  this->make_copy(another);
}

//Destructor for the RcppParallelBatchSimulator class.
RcppParallelBatchSimulator::~RcppParallelBatchSimulator(void)
{
}

void RcppParallelBatchSimulator::operator=(const RcppParallelBatchSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }

  BatchSimulator::operator=(another);
  this->make_copy(another);
}

BatchSimulator* RcppParallelBatchSimulator::duplicate(void)const
{
  return( new RcppParallelBatchSimulator(*this));
}

void RcppParallelBatchSimulator::make_copy(const RcppParallelBatchSimulator &another)
{

}

Particles RcppParallelBatchSimulator::simulate(void) const
{
  return Particles();
}
