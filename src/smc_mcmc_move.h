#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

#include "smc.h"

class SMCMCMCMove : public SMC
{
public:

  SMCMCMCMove();

  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              const Data* data_in,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in);
  SMCMCMCMove(const SMCMCMCMove &another);
  virtual ~SMCMCMCMove(void);

  void operator=(const SMCMCMCMove &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;

protected:

  void smc_update(SMCOutput* current_state);

  void make_copy(const SMCMCMCMove &another);

  //void smc_step();

  //void weight_update();
};

#endif
