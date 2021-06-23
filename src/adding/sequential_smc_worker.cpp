#include "sequential_smc_worker.h"

//Default constructor.
SequentialSMCWorker::SequentialSMCWorker(void)
  :SMCWorker()
{
}

//Copy constructor for the SequentialSMCWorker class.
SequentialSMCWorker::SequentialSMCWorker(const SequentialSMCWorker &another)
  :SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialSMCWorker class.
SequentialSMCWorker::~SequentialSMCWorker(void)
{
}

void SequentialSMCWorker::operator=(const SequentialSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* SequentialSMCWorker::duplicate(void)const
{
  return( new SequentialSMCWorker(*this));
}

void SequentialSMCWorker::make_copy(const SequentialSMCWorker &another)
{

}

Particles SequentialSMCWorker::simulate(void) const
{
  return Particles();
}
