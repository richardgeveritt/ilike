#include "model_and_algorithm.h"

#include <array>

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
  this->is_simulate_methods = another.is_simulate_methods;
  this->simulate_priors = another.simulate_priors;
}
