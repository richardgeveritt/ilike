#include "sequential_smc_worker.h"
#include "particle_simulator.h"
#include "smc.h"

//Default constructor.
SequentialSMCWorker::SequentialSMCWorker(void)
  :SMCWorker()
{
}

SequentialSMCWorker::SequentialSMCWorker(SMC* the_smc_in,
                                         ParticleSimulator* particle_simulator_in)
  :SMCWorker(the_smc_in, particle_simulator_in)
{
  this->simulate_output = std::vector<Particle>(this->get_number_of_particles());
}

//Copy constructor for the SequentialSMCWorker class.
SequentialSMCWorker::SequentialSMCWorker(const SequentialSMCWorker &another)
  :SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialSMCWorker class.
SequentialSMCWorker::~SequentialSMCWorker(void)
{
}

void SequentialSMCWorker::operator=(const SequentialSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* SequentialSMCWorker::duplicate(void)const
{
  return( new SequentialSMCWorker(*this));
}

void SequentialSMCWorker::make_copy(const SequentialSMCWorker &another)
{

}

std::vector<Particle> SequentialSMCWorker::get_particles() const
{
  return this->simulate_output;
}

void SequentialSMCWorker::specific_simulate()
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->simulate_output[i] = (*this->particle_simulator)(local_rng);
  }
}

void SequentialSMCWorker::simulate_and_weight()
{
}
