#include "model_and_algorithm.h"
#include "likelihood_estimator.h"
#include "importance_sampler.h"
#include "parameter_particle_simulator.h"
#include "exact_likelihood_estimator.h"

//using namespace std::placeholders;

// Particle simulate_particle_from_simulate_parameters(SimulateDistributionPtr simulate_distribution, RandomNumberGenerator &rng)
// {
//   return Particle(simulate_distribution(rng));
// }

ModelAndAlgorithm::ModelAndAlgorithm()
{
  this->likelihood_estimators.resize(0);
  this->particle_simulator = NULL;
  this->observed_data = NULL;
}

ModelAndAlgorithm::ModelAndAlgorithm(const ModelAndAlgorithm &another)
{
  this->make_copy(another);
}

ModelAndAlgorithm::~ModelAndAlgorithm(void)
{
  for (std::vector<LikelihoodEstimator*>::iterator i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  if (this->particle_simulator!=NULL)
    delete this->particle_simulator;
}

void ModelAndAlgorithm::operator=(const ModelAndAlgorithm &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void ModelAndAlgorithm::make_copy(const ModelAndAlgorithm &another)
{
  //this->is_simulate_methods = another.is_simulate_methods;
  //this->simulate_priors = another.simulate_priors;
  this->likelihood_estimators.resize(0);
  this->likelihood_estimators.reserve(another.likelihood_estimators.size());
  for (std::vector<LikelihoodEstimator*>::const_iterator i=another.likelihood_estimators.begin();
       i!=another.likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimators.push_back((*i)->duplicate());
  }
  if (another.particle_simulator!=NULL)
    this->particle_simulator = another.particle_simulator->duplicate();
  this->observed_data = another.observed_data;
  //this->initial_seed = another.initial_seed;
}

// ImportanceSampler* ModelAndAlgorithm::get_is(RandomNumberGenerator* rng_in,
//                                              size_t* seed_in,
//                                              bool parallel_in,
//                                              const Data* data_in,
//                                              size_t number_of_particles_in,
//                                              SimulateDistributionPtr simulate_distribution_in,
//                                              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in)
// {
//
//
//
//
//   //auto thing = std::bind(simulate_particle_from_simulate_parameters, simulate_distribution_in, _2);
//
//   return new ImportanceSampler(rng_in, seed_in, parallel_in, data_in, number_of_particles_in, *this);
// }
