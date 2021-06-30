#ifndef RCPPPARALLELSMCWORKER_H
#define RCPPPARALLELSMCWORKER_H

#include <RcppParallel.h>

#include "smc_worker.h"
#include "distributions.h"

class RcppParallelSMCWorker;
class ParticleSimulator;

class SimulateWorker : public RcppParallel::Worker {

public:

  SimulateWorker();

  SimulateWorker(RcppParallelSMCWorker* smc_worker_in,
                 const ParticleSimulator* particle_simulator_in,
                 const size_t num_particles_in);

  ~SimulateWorker();

  SimulateWorker(const SimulateWorker &another);

  void operator=(const SimulateWorker &another);

  void operator()(std::size_t begin, std::size_t end);

  std::vector<Particle> simulate_output;

private:

  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;

  const ParticleSimulator* particle_simulator;

  void make_copy(const SimulateWorker &another);

  //uint64_t seed;

  //std::vector<Parameters> input;

  //std::vector<double> dummy_input;

  //EvaluateLogDistributionPtr evaluate_log_prior;

  //SimulateDistributionPtr simulate_prior;

  //random_number_generator rng;

};


class RcppParallelSMCWorker : public SMCWorker
{
public:

  RcppParallelSMCWorker(uint64_t seed_in,
                        const ParticleSimulator* particle_simulator_in,
                        const size_t num_particles_in);
  virtual ~RcppParallelSMCWorker();

  RcppParallelSMCWorker(const RcppParallelSMCWorker &another);
  void operator=(const RcppParallelSMCWorker &another);
  SMCWorker* duplicate() const;

  Particles simulate() const;

protected:

  void make_copy(const RcppParallelSMCWorker &another);

  uint64_t seed;
  RandomNumberGenerator rng;

  friend class SimulateWorker;
  SimulateWorker simulate_worker;

};

#endif
