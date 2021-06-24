#include "smc.h"

#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

class ImportanceSampler : public SMC
{
public:

  ImportanceSampler(const ModelAndAlgorithm* model_and_algorithm_in);
  ImportanceSampler(const ImportanceSampler &another);
  virtual ~ImportanceSampler(void);

  void operator=(const ImportanceSampler &another);
  SMC* duplicate() const;

  SMCOutput do_smc();

protected:

  void make_copy(const ImportanceSampler &another);

  void smc_step();

  void weight_update();

  //double single_particle_weight_update() const;
};

#endif
