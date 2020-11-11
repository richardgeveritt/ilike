#include "smc.h"

#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

class SMCMCMCMove : public SMC
{
public:

  SMCMCMCMove(void);
  SMCMCMCMove(const SMCMCMCMove &another);
  virtual ~SMCMCMCMove(void);

  void operator=(const SMCMCMCMove &another);
  SMC* duplicate() const;

protected:

  void make_copy(const SMCMCMCMove &another);
};

#endif
