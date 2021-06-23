#include "smc_worker.h"

#include <vector>

//#include "simulator.h"

#ifndef SEQUENTIALSMCWORKER_H
#define SEQUENTIALSMCWORKER_H

class SequentialSMCWorker : public SMCWorker
{
public:

  SequentialSMCWorker(void);
  SequentialSMCWorker(const SequentialSMCWorker &another);
  virtual ~SequentialSMCWorker(void);

  void operator=(const SequentialSMCWorker &another);
  SMCWorker* duplicate() const;

  Particles simulate(void) const;

protected:

  void make_copy(const SequentialSMCWorker &another);

  // Don't think we want to store this here, since found in particles (I think).
  //std::vector<Simulator*> simulate_priors;
  //std::vector<> simulate_for_likelihoods;

};

#endif
