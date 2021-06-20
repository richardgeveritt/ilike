#include "batch_simulator.h"

#include <vector>

//#include "simulator.h"

#ifndef SEQUENTIALBATCHSIMULATOR_H
#define SEQUENTIALBATCHSIMULATOR_H

class SequentialBatchSimulator : public BatchSimulator
{
public:

  SequentialBatchSimulator(void);
  SequentialBatchSimulator(const SequentialBatchSimulator &another);
  virtual ~SequentialBatchSimulator(void);

  void operator=(const SequentialBatchSimulator &another);
  BatchSimulator* duplicate() const;

  Particles simulate(void) const;

protected:

  void make_copy(const SequentialBatchSimulator &another);

  // Don't think we want to store this here, since found in particles (I think).
  //std::vector<Simulator*> simulate_priors;
  //std::vector<> simulate_for_likelihoods;

};

#endif
