#include "likelihood_estimator.h"
#include "model_and_algorithm.h"
#include "likelihood_estimator_output.h"
#include "smc_worker.h"

LikelihoodEstimator::LikelihoodEstimator()
{
}

LikelihoodEstimator::LikelihoodEstimator(RandomNumberGenerator* rng_in,
                                         size_t* seed_in,
                                         const Data* data_in)
{
  this->data = data_in;
  this->model_and_algorithm = ModelAndAlgorithm() ;
  this->rng = rng_in;
  this->seed = seed_in;
}

LikelihoodEstimator::~LikelihoodEstimator()
{

}

LikelihoodEstimator::LikelihoodEstimator(const LikelihoodEstimator &another)
{
  this->make_copy(another);
}

void LikelihoodEstimator::operator=(const LikelihoodEstimator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void LikelihoodEstimator::make_copy(const LikelihoodEstimator &another)
{
  this->data = another.data;
  this->model_and_algorithm = another.model_and_algorithm;
  this->rng = another.rng;
  this->seed = another.seed;
}
