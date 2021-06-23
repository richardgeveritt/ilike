#include "smc_worker.h"

SMCWorker::SMCWorker(void)
{
}

SMCWorker::SMCWorker(const SMCWorker &another)
{
  this->make_copy(another);
}

SMCWorker::~SMCWorker(void)
{
}

void SMCWorker::operator=(const SMCWorker &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMCWorker::make_copy(const SMCWorker &another)
{
  //Does nothing since no member variables to copy.
}
