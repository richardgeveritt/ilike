#ifndef MODELANDALGORITHM_H
#define MODELANDALGORITHM_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "distributions.h"
#include "function_pointers.h"


class LikelihoodEstimator;
class ImportanceSampler;
class ParticleSimulator;

//typedef Particle (*SimulateParticlePtr)(random_number_generator &rng);



class ModelAndAlgorithm
{

public:

  //std::vector<std::string> is_simulate_methods;

  //BatchSimulator* simulate_priors;

  //Particle simulate(random_number_generator &rng) const;
  //std::vector<> simulate_for_likelihoods;

  ModelAndAlgorithm();

  virtual ~ModelAndAlgorithm();

  ModelAndAlgorithm(const ModelAndAlgorithm &another);
  void operator=(const ModelAndAlgorithm &another);
//
//   ImportanceSampler* get_is(RandomNumberGenerator* rng_in,
//                             size_t* seed_in,
//                             bool parallel_in,
//                             const Data* data_in,
//                             size_t number_of_particles_in,
//                             SimulateDistributionPtr simulate_distribution_in,
//                             EvaluateLogLikelihoodPtr evaluate_log_likelihood_in);

  // Stored here.
  std::vector<LikelihoodEstimator*> likelihood_estimators;

  // Stored here.
  ParticleSimulator* particle_simulator;

  // Not stored here. Stored in "main".
  const Data* observed_data;

protected:

  void make_copy(const ModelAndAlgorithm &another);



  //size_t initial_seed;

};

#endif
