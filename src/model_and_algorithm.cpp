#include "model_and_algorithm.h"
#include "likelihood_estimator.h"
#include "importance_sampler.h"
//#include "parameter_particle_simulator.h"

//using namespace std::placeholders;

// Particle simulate_particle_from_simulate_parameters(SimulateDistributionPtr simulate_distribution, RandomNumberGenerator &rng)
// {
//   return Particle(simulate_distribution(rng));
// }

ModelAndAlgorithm::ModelAndAlgorithm()
{
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
    this->likelihood_estimators.push_back((*i)->duplicate());
  }
}

ImportanceSampler* ModelAndAlgorithm::GetIS(const SimulateDistributionPtr simulate_distribution_in,
                              const EvaluateLogLikelihoodPtr evaluate_log_likelihood_in)
{
  // Need to construct LikelihoodEstimator to read in to this constructor.
  //ParameterParticleSimulator parameter_particle_simulator(simulate_distribution_in);

  //auto thing = std::bind(simulate_particle_from_simulate_parameters, simulate_distribution_in, _2);

  return new ImportanceSampler();
}
