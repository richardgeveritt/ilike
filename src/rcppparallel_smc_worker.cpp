#include "rcppparallel_smc_worker.h"
#include "particle_simulator.h"
#include "smc.h"

SimulateWorker::SimulateWorker()
{
}

SimulateWorker::SimulateWorker(RcppParallelSMCWorker* smc_worker_in)
{
  this->smc_worker = smc_worker_in;
  this->simulate_output = std::vector<Particle>(this->smc_worker->get_number_of_particles());
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
  RandomNumberGenerator local_rng(*this->smc_worker->get_rng());
  local_rng.seed(this->smc_worker->get_seed(),end);
  for (std::size_t i = begin; i < end; ++i)
  {
    this->simulate_output[i] = (*this->smc_worker->particle_simulator)(local_rng);
  }
}

void SimulateWorker::make_copy(const SimulateWorker &another)
{
  this->smc_worker = another.smc_worker;
  this->simulate_output = another.simulate_output;
}

//Default constructor.
RcppParallelSMCWorker::RcppParallelSMCWorker()
{
}

RcppParallelSMCWorker::RcppParallelSMCWorker(SMC* the_smc_in,
                                             ParticleSimulator* particle_simulator_in,
                                             size_t grain_size_in)
  :SMCWorker(the_smc_in, particle_simulator_in)
{
  this->simulate_worker = SimulateWorker(this);
  this->grain_size = grain_size_in;
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
  //this->seed = another.seed;
  //this->rng = another.rng;
  this->simulate_worker = another.simulate_worker;
  this->grain_size = another.grain_size;
  //this->particle_simulator = another.particle_simulator->duplicate();
}

std::vector<Particle> RcppParallelSMCWorker::get_particles() const
{
  return this->simulate_worker.simulate_output;
}

void RcppParallelSMCWorker::specific_simulate()
{
  parallelFor(0, this->get_number_of_particles(), this->simulate_worker, this->grain_size);
}

void RcppParallelSMCWorker::simulate_and_weight()
{

}
