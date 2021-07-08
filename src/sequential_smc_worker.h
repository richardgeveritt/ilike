#ifndef SEQUENTIALSMCWORKER_H
#define SEQUENTIALSMCWORKER_H

#include "smc_worker.h"

#include <vector>

#include "particles.h"

class ParticleSimulator;

class SequentialSMCWorker : public SMCWorker
{
public:

  SequentialSMCWorker(void);
  virtual ~SequentialSMCWorker(void);

  SequentialSMCWorker(SMC* the_smc_in,
                      ParticleSimulator* particle_simulator_in);

  SequentialSMCWorker(const SequentialSMCWorker &another);
  void operator=(const SequentialSMCWorker &another);
  SMCWorker* duplicate() const;

  void simulate_and_weight(void);

  std::vector<Particle> get_particles() const;

protected:

  void specific_simulate();

  void make_copy(const SequentialSMCWorker &another);

  // Don't think we want to store this here, since found in particles (I think).
  //std::vector<Simulator*> simulate_priors;
  //std::vector<> simulate_for_likelihoods;

  std::vector<Particle> simulate_output;

};

#endif
