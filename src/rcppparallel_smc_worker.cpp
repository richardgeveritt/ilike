#include "rcppparallel_smc_worker.h"
#include "particle_simulator.h"

SimulateWorker::SimulateWorker()
{
}

SimulateWorker::SimulateWorker(RcppParallelSMCWorker* smc_worker_in,
               const ParticleSimulator* particle_simulator_in,
               const size_t num_particles_in)
{
  this->smc_worker = smc_worker_in;
  this->particle_simulator = particle_simulator_in;
  this->simulate_output = std::vector<Particle>(num_particles_in);
}

SimulateWorker::~SimulateWorker()
{

}

SimulateWorker::SimulateWorker(const SimulateWorker &another)
{
  this->make_copy(another);
}

void SimulateWorker::operator=(const SimulateWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->make_copy(another);
}

void SimulateWorker::operator()(std::size_t begin, std::size_t end)
{
  RandomNumberGenerator local_rng(this->smc_worker->rng);
  local_rng.seed(this->smc_worker->seed,end);
  for (std::size_t i = begin; i < end; ++i)
  {
    this->simulate_output[i] = (*this->particle_simulator)(this->smc_worker->rng);
  }
}

void SimulateWorker::make_copy(const SimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->particle_simulator = another.particle_simulator;//->duplicate();
}

//Default constructor.
RcppParallelSMCWorker::RcppParallelSMCWorker(uint64_t seed_in,
                                             const ParticleSimulator* particle_simulator_in,
                                             const size_t num_particles_in)
  :SMCWorker(), simulate_worker()
{
  this->seed = seed_in;
  this->simulate_worker = SimulateWorker(this,
                                         particle_simulator_in,
                                         num_particles_in);
}

//Copy constructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::RcppParallelSMCWorker(const RcppParallelSMCWorker &another)
  :SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the RcppParallelSMCWorker class.
RcppParallelSMCWorker::~RcppParallelSMCWorker(void)
{
  // if (this->particle_simulator!=NULL)
  //   delete this->particle_simulator;
}

void RcppParallelSMCWorker::operator=(const RcppParallelSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* RcppParallelSMCWorker::duplicate(void)const
{
  return( new RcppParallelSMCWorker(*this));
}

void RcppParallelSMCWorker::make_copy(const RcppParallelSMCWorker &another)
{
  this->seed = another.seed;
  this->rng = another.rng;
  this->simulate_worker = another.simulate_worker;
}

Particles RcppParallelSMCWorker::simulate(void) const
{
  return Particles();
}
