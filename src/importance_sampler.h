#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "smc.h"
#include "model_and_algorithm.h"
#include "distributions.h"
#include "parameters.h"

class ImportanceSampler : public SMC
{
public:

  ImportanceSampler();
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    bool parallel_in,
                    const Data* data_in,
                    size_t number_of_particles_in,
                    SimulateDistributionPtr simulate_distribution_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in);

  ImportanceSampler(const ImportanceSampler &another);
  virtual ~ImportanceSampler(void);

  void operator=(const ImportanceSampler &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;

  LikelihoodEstimatorOutput* run();
  LikelihoodEstimatorOutput* run(const Parameters &parameters);

protected:

  void smc_update(SMCOutput* current_state);

  void make_copy(const ImportanceSampler &another);

  // void smc_step();
  //
  // void weight_update();

  //double single_particle_weight_update() const;
};

#endif
