#include "sequencer.h"

#ifndef DETERMINISTICSEQUENCER_H
#define DETERMINISTICSEQUENCER_H

class DeterministicSequencer : public Sequencer
{
public:

  DeterministicSequencer(void);
  DeterministicSequencer(const DeterministicSequencer &another);
  virtual ~DeterministicSequencer(void);

  void operator=(const DeterministicSequencer &another);
  Sequencer* duplicate() const;

protected:

  void make_copy(const DeterministicSequencer &another);
};

#endif
