#include "smc_worker.h"
#include "smc.h"
#include "particle_simulator.h"

SMCWorker::SMCWorker(void)
{
  this->particle_simulator = NULL;
}

SMCWorker::SMCWorker(SMC* the_smc_in,
                     ParticleSimulator* particle_simulator_in)
{
  this->the_smc = the_smc_in;
  this->particle_simulator = particle_simulator_in;
}

SMCWorker::SMCWorker(const SMCWorker &another)
{
  this->make_copy(another);
}

SMCWorker::~SMCWorker(void)
{
}

void SMCWorker::operator=(const SMCWorker &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMCWorker::make_copy(const SMCWorker &another)
{
  this->the_smc = another.the_smc;
}

void SMCWorker::simulate()
{
  this->specific_simulate();
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}

size_t SMCWorker::get_number_of_particles() const
{
  return this->the_smc->number_of_particles;
}

RandomNumberGenerator* SMCWorker::get_rng()
{
  return this->the_smc->rng;
}

size_t SMCWorker::get_seed() const
{
  return *this->the_smc->seed;
}

void SMCWorker::set_seed(size_t seed_in)
{
  *this->the_smc->seed = seed_in;
}
