#include "smc_worker.h"
#include "smc.h"
#include "particle_simulator.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

namespace ilike
{
SMCWorker::SMCWorker()
{
}

// Constructor for IS where prior is proposal.
SMCWorker::SMCWorker(SMC* the_smc_in)
{
  this->the_smc = the_smc_in;
}

SMCWorker::SMCWorker(const SMCWorker &another)
{
  this->make_copy(another);
}

SMCWorker::~SMCWorker()
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

void SMCWorker::simulate(Particles* next_particles)
{
  this->specific_simulate(next_particles);
  
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
}

void SMCWorker::simulate(Particles* next_particles,
                         const Parameters &conditioned_on_parameters)
{
  this->specific_simulate(next_particles,
                          conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
}

void SMCWorker::move(Particles* next_particles,
                     const Particles* current_particles)
{
  this->specific_move(next_particles,
                      current_particles);
  
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}

void SMCWorker::subsample_simulate(Particles* next_particles,
                                   const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_simulate(next_particles,
                                    conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
}

void SMCWorker::subsample_move(Particles* next_particles,
                               const Particles* current_particles)
{
  this->subsample_specific_move(next_particles,
                                current_particles);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
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
}
