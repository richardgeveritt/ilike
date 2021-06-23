#include "rcppparallel_smc_worker.h"

//Default constructor.
RcppParallelSMCWorker::RcppParallelSMCWorker(void)
  :SMCWorker()
{
}

//Copy constructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::RcppParallelSMCWorker(const RcppParallelSMCWorker &another)
  :SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::~RcppParallelSMCWorker(void)
{
}

void RcppParallelSMCWorker::operator=(const RcppParallelSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* RcppParallelSMCWorker::duplicate(void)const
{
  return( new RcppParallelSMCWorker(*this));
}

void RcppParallelSMCWorker::make_copy(const RcppParallelSMCWorker &another)
{

}

Particles RcppParallelSMCWorker::simulate(void) const
{
  return Particles();
}
