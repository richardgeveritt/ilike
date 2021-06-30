#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

#include "smc.h"

class SMCMCMCMove : public SMC
{
public:

  SMCMCMCMove(const ModelAndAlgorithm &model_and_algorithm_in,
              const Data* data_in);
  SMCMCMCMove(const SMCMCMCMove &another);
  virtual ~SMCMCMCMove(void);

  void operator=(const SMCMCMCMove &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;

  SMCOutput do_smc();

protected:

  void make_copy(const SMCMCMCMove &another);

  void smc_step();

  void weight_update();
};

#endif
