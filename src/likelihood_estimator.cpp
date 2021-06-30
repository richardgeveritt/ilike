#include "likelihood_estimator.h"
#include "model_and_algorithm.h"

LikelihoodEstimator::LikelihoodEstimator()
{
}

LikelihoodEstimator::LikelihoodEstimator(const ModelAndAlgorithm &model_and_algorithm_in,
                                         const Data* observed_data_in)
{
  this->observed_data = observed_data_in;
  this->model_and_algorithm = model_and_algorithm_in;
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
  this->observed_data = another.observed_data;
  this->model_and_algorithm = another.model_and_algorithm;
}
