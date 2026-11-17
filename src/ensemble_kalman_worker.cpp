#include "ensemble_kalman_worker.h"
#include "ensemble_kalman.h"
#include "particle_simulator.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

namespace ilike
{
EnsembleKalmanWorker::EnsembleKalmanWorker()
{
}

// Constructor for IS where prior is proposal.
EnsembleKalmanWorker::EnsembleKalmanWorker(EnsembleKalman* the_enk_in)
{
  this->the_enk = the_enk_in;
}

EnsembleKalmanWorker::EnsembleKalmanWorker(const EnsembleKalmanWorker &another)
{
  this->make_copy(another);
}

EnsembleKalmanWorker::~EnsembleKalmanWorker()
{
}

void EnsembleKalmanWorker::operator=(const EnsembleKalmanWorker &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void EnsembleKalmanWorker::make_copy(const EnsembleKalmanWorker &another)
{
  this->the_enk = another.the_enk;

}

void EnsembleKalmanWorker::simulate(Ensemble* next_ensemble,
                                    const Index* index)
{
  this->specific_simulate(next_ensemble,
                          index);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  //this->the_enk->rng->seed(this->get_seed());
  //this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

void EnsembleKalmanWorker::simulate(Ensemble* next_ensemble,
                                    const Index* index,
                                    const Parameters &conditioned_on_parameters)
{
  this->specific_simulate(next_ensemble,
                          index,
                          conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  //this->the_enk->rng->seed(this->get_seed());
  //this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
  
}

void EnsembleKalmanWorker::move(Ensemble* next_particles,
                                Ensemble* current_particles)
{
  this->specific_move(next_particles,
                      current_particles);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

void EnsembleKalmanWorker::subsample_move(Ensemble* next_particles,
                                          Ensemble* current_particles)
{
  this->subsample_specific_move(next_particles,
                                current_particles);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

size_t EnsembleKalmanWorker::get_number_of_ensemble_members() const
{
  return this->the_enk->number_of_ensemble_members;
}

RandomNumberGenerator* EnsembleKalmanWorker::get_rng()
{
  return this->the_enk->rng;
}

size_t EnsembleKalmanWorker::get_seed() const
{
  return *this->the_enk->seed;
}

void EnsembleKalmanWorker::set_seed(size_t seed_in)
{
  *this->the_enk->seed = seed_in;
}

void EnsembleKalmanWorker::set_enk(EnsembleKalman* the_enk_in)
{
  this->the_enk = the_enk_in;
}
}
