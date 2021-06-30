#ifndef MODELANDALGORITHM_H
#define MODELANDALGORITHM_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "simulation.h"
#include "function_pointers.h"


class LikelihoodEstimator;


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

  void SetIS(const SimulateDistributionPtr simulate_distribution_in,
             const EvaluateLogLikelihoodPtr evaluate_log_likelihood_in);

protected:

  void make_copy(const ModelAndAlgorithm &another);

  // Stored here.
  std::vector<LikelihoodEstimator*> likelihood_estimators;

};

#endif
