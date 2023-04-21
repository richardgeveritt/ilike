#include "smc_worker.h"
#include "smc.h"
#include "particle_simulator.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

SMCWorker::SMCWorker(void)
{
  //this->particle_simulator = NULL;
}

// Constructor for IS where prior is proposal.
SMCWorker::SMCWorker(SMC* the_smc_in)
{
  this->the_smc = the_smc_in;
  //this->particle_simulator = particle_simulator_in;
  //this->factors = factors_in;
}

SMCWorker::SMCWorker(const SMCWorker &another)
{
  this->make_copy(another);
}

SMCWorker::~SMCWorker()
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

void SMCWorker::operator=(const SMCWorker &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMCWorker::make_copy(const SMCWorker &another)
{
  this->the_smc = another.the_smc;
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

void SMCWorker::simulate(Particles* next_particles)
{
  this->specific_simulate(next_particles);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
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
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
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
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}

/*
void SMCWorker::move(Particles* next_particles,
                     const Particles* current_particles,
                     const Parameters &conditioned_on_parameters)
{
  this->specific_move(next_particles,
                      current_particles,
                      conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}
*/

void SMCWorker::subsample_simulate(Particles* next_particles,
                                   const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_simulate(next_particles,
                          conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Initialise weight.
  //this->simulated_particles().initialise_weights();
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
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
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}

/*
void SMCWorker::subsample_move(Particles* next_particles,
                               const Particles* current_particles,
                               const Parameters &conditioned_on_parameters)
{
  this->subsample_specific_move(next_particles,
                                current_particles,
                                conditioned_on_parameters);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
  
  // Simulate random numbers for resampling.
  this->the_smc->rng->seed(this->get_seed(),this->get_number_of_particles());
  next_particles->simulate_resampling_variables(*this->the_smc->rng);
  this->set_seed(this->get_seed() + this->get_number_of_particles());
}
*/

//void SMCWorker::weight(const Parameters &conditioned_on_parameters)
//{
//  this->specific_weight(conditioned_on_parameters);
//}

//void SMCWorker::simulate_and_weight()
//{
//  this->specific_simulate_and_weight();
//  this->set_seed(this->get_seed() + this->get_number_of_particles());
//}

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
