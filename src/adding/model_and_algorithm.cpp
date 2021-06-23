#include "model_and_algorithm.h"

#include <array>

using namespace std::placeholders;

Particle simulate_particle_from_simulate_parameters(SimulateDistributionPtr simulate_distribution, random_number_generator &rng)
{
  return Particle(simulate_distribution(rng));
}

ModelAndAlgorithm::ModelAndAlgorithm()
{
}

ModelAndAlgorithm::ModelAndAlgorithm(const ModelAndAlgorithm &another)
{
  this->make_copy(another);
}

ModelAndAlgorithm::~ModelAndAlgorithm(void)
{
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
}

void ModelAndAlgorithm::SetIS(const SimulateDistributionPtr simulate_distribution_in,
                              const EvaluateLogLikelihoodPtr evaluate_log_likelihood_in)
{
  auto thing = std::bind(simulate_particle_from_simulate_parameters, simulate_distribution_in, _2);
}
