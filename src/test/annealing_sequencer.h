#include "sequencer.h"

#ifndef ANNEALINGSEQUENCER_H
#define ANNEALINGSEQUENCER_H

class AnnealingSequencer : public Sequencer
{
public:

  AnnealingSequencer(void);
  AnnealingSequencer(const AnnealingSequencer &another);
  virtual ~AnnealingSequencer(void);

  void operator=(const AnnealingSequencer &another);
  Sequencer* duplicate() const;

protected:

  void make_copy(const AnnealingSequencer &another);
};

#endif
