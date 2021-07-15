#ifndef RCPPPARALLELSMCWORKER_H
#define RCPPPARALLELSMCWORKER_H

#include <RcppParallel.h>

#include "smc_worker.h"
#include "distributions.h"

class SMC;
class RcppParallelSMCWorker;
class ParticleSimulator;

class SimulateWorker : public RcppParallel::Worker {

public:

  SimulateWorker();

  SimulateWorker(RcppParallelSMCWorker* smc_worker_in);

  ~SimulateWorker();

  SimulateWorker(const SimulateWorker &another);

  void operator=(const SimulateWorker &another);

  void operator()(std::size_t begin, std::size_t end);

  std::vector<Particle> simulate_output;

private:

  friend class RcppParallelSMCWorker;
  RcppParallelSMCWorker* smc_worker;

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

  RcppParallelSMCWorker();
  RcppParallelSMCWorker(SMC* the_smc,
                        ParticleSimulator* particle_simulator_in,
                        size_t grain_size_in);
  virtual ~RcppParallelSMCWorker();

  RcppParallelSMCWorker(const RcppParallelSMCWorker &another);
  void operator=(const RcppParallelSMCWorker &another);
  SMCWorker* duplicate() const;

  std::vector<Particle> get_particles() const;

  // Simulate from the proposal and weight.
  void simulate_and_weight();

protected:

  // Just does the simulation from the proposal.
  void specific_simulate();

  void make_copy(const RcppParallelSMCWorker &another);

  //uint64_t seed;
  //RandomNumberGenerator rng;

  size_t grain_size;

  friend class SimulateWorker;
  SimulateWorker simulate_worker;

};

#endif
