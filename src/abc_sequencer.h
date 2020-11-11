#include "sequencer.h"

#ifndef ABCSEQUENCER_H
#define ABCSEQUENCER_H

class ABCSequencer : public Sequencer
{
public:

  ABCSequencer(void);
  ABCSequencer(const ABCSequencer &another);
  virtual ~ABCSequencer(void);

  void operator=(const ABCSequencer &another);
  Sequencer* duplicate() const;

protected:

  void make_copy(const ABCSequencer &another);
};

#endif
