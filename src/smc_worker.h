#ifndef SMCWORKER_H
#define SMCWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "particles.h"
#include "distributions.h"

class SMC;
class ParticleSimulator;

class SMCWorker
{
public:

  SMCWorker();
  SMCWorker(SMC* the_smc_in,
            ParticleSimulator* particle_simulator_in);
  virtual ~SMCWorker(void);

  SMCWorker(const SMCWorker &another);
  void operator=(const SMCWorker &another);
  virtual SMCWorker* duplicate() const=0;

  size_t get_number_of_particles() const;
  RandomNumberGenerator* get_rng();
  size_t get_seed() const;
  void set_seed(size_t seed_in);

  void simulate();

  virtual void simulate_and_weight()=0;

  virtual std::vector<Particle> get_particles() const=0;

protected:

  virtual void specific_simulate()=0;

  void make_copy(const SMCWorker &another);

  SMC* the_smc; // not stored here

  // stored in model_and_algorithm
  ParticleSimulator* particle_simulator;

  // Add in function to do simuulation of parts that cannot be done in parallel.

};

#endif
