#include "batch_simulator.h"

#ifndef RCPPPARALLELBATCHSIMULATOR_H
#define RCPPPARALLELBATCHSIMULATOR_H

class RcppParallelBatchSimulator : public BatchSimulator
{
public:

  RcppParallelBatchSimulator(void);
  RcppParallelBatchSimulator(const RcppParallelBatchSimulator &another);
  virtual ~RcppParallelBatchSimulator(void);

  void operator=(const RcppParallelBatchSimulator &another);
  BatchSimulator* duplicate() const;

  Particles simulate(void) const;

protected:

  void make_copy(const RcppParallelBatchSimulator &another);

};

#endif
