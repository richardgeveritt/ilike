#include "smc_worker.h"

#ifndef RCPPPARALLELSMCWORKER_H
#define RCPPPARALLELSMCWORKER_H

class RcppParallelSMCWorker : public SMCWorker
{
public:

  RcppParallelSMCWorker();
  RcppParallelSMCWorker(const RcppParallelSMCWorker &another);
  virtual ~RcppParallelSMCWorker();

  void operator=(const RcppParallelSMCWorker &another);
  SMCWorker* duplicate() const;

  Particles simulate() const;

protected:

  void make_copy(const RcppParallelSMCWorker &another);

  uint64_t seed;
  random_number_generator rng;

};

#endif
