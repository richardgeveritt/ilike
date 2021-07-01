#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "smc.h"

class ImportanceSampler : public SMC
{
public:

  ImportanceSampler();
  ImportanceSampler(const ModelAndAlgorithm &model_and_algorithm_in,
                    const Data* data_in);
  ImportanceSampler(const ImportanceSampler &another);
  virtual ~ImportanceSampler(void);

  void operator=(const ImportanceSampler &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;

protected:

  void make_copy(const ImportanceSampler &another);

  void smc_step();

  void weight_update();

  //double single_particle_weight_update() const;
};

#endif
