#include "ensemble_kalman_worker.h"
#include "ensemble_kalman.h"
#include "particle_simulator.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

EnsembleKalmanWorker::EnsembleKalmanWorker()
{
  //this->particle_simulator = NULL;
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
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::iterator i=this->likelihood_estimator_outputs.begin();
    //   i!=this->likelihood_estimator_outputs.end();
    //   ++i)
  //{
  //  for (std::vector<LikelihoodEstimatorOutput*>::iterator j=i->begin();
  //       j!=i->end();
  //       ++j)
  //  {
  //    if (*j!=NULL)
  //      delete *j;
  //  }
  //}
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
  //this->particle_simulator = another.particle_simulator;
  //this->likelihood_estimators = another.likelihood_estimators;
  //this->output = another.output;
  
  //this->likelihood_estimator_outputs.resize(0);
  //this->likelihood_estimator_outputs.reserve(another.likelihood_estimator_outputs.size());
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator i=this->likelihood_estimator_outputs.begin();
  //     i!=this->likelihood_estimator_outputs.end();
  //     ++i)
  //{
  //  std::vector<LikelihoodEstimatorOutput*> inner_vector;
  //  inner_vector.reserve(this->get_number_of_particles());
    
  //  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator j=i->begin();
  //       j!=i->end();
  //       ++j)
  //  {
  //    if (*j!=NULL)
  //      inner_vector.push_back((*j)->duplicate());
  //    else
  //      inner_vector.push_back(NULL);
  //  }
  //  this->likelihood_estimator_outputs.push_back(inner_vector);
  //}
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
  //this->the_enk->rng->seed(this->get_seed(),this->get_number_of_ensemble_members());
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
  //this->the_enk->rng->seed(this->get_seed(),this->get_number_of_ensemble_members());
  //this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
  
}

void EnsembleKalmanWorker::move(Ensemble* next_particles,
                                Ensemble* current_particles)
{
  this->specific_move(next_particles,
                      current_particles);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

void EnsembleKalmanWorker::move(Ensemble* next_particles,
                                Ensemble* current_particles,
                                const Parameters &conditioned_on_parameters)
{
  this->specific_move(next_particles,
                      current_particles,
                      conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

void EnsembleKalmanWorker::subsample_move(Ensemble* next_particles,
                                          Ensemble* current_particles,
                                          const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_move(next_particles,
                                current_particles,
                                conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_ensemble_members());
}

/*
void EnsembleKalmanWorker::subsample_simulate(Particles* next_particles,
                                   const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_simulate(next_particles,
                          conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_enk->rng->seed(this->get_seed(),this->get_number_of_particles());
  next_particles->simulate_resampling_variables(*this->the_enk->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
}

void EnsembleKalmanWorker::subsample_move(Particles* next_particles,
                               const Particles* current_particles,
                               const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_move(next_particles,
                                current_particles,
                                conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}
*/

//void EnsembleKalmanWorker::weight(const Parameters &conditioned_on_parameters)
//{
//  this->specific_weight(conditioned_on_parameters);
//}

//void EnsembleKalmanWorker::simulate_and_weight()
//{
//  this->specific_simulate_and_weight();
//  this->set_seed(this->get_seed() + this->get_number_of_particles());
//}

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
