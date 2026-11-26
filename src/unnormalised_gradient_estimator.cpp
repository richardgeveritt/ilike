#include "unnormalised_gradient_estimator.h"
#include "unnormalised_gradient_estimator_output.h"

namespace ilike
{

UnnormalisedGradientEstimator::UnnormalisedGradientEstimator()
  :LikelihoodEstimator()
{
  this->statistic = NULL;
  this->exact = NULL;
  this->sweep = NULL;
  this->sweeps = 0;
  this->points = 1;
}

UnnormalisedGradientEstimator::UnnormalisedGradientEstimator(RandomNumberGenerator* rng_in,
                                                             size_t* seed_in,
                                                             Data* data_in,
                                                             const std::string &variable_in,
                                                             UnnormalisedPotentialPtr statistic_in,
                                                             UnnormalisedExactFunctionPtr exact_in,
                                                             UnnormalisedSweepFunctionPtr sweep_in,
                                                             size_t sweeps_in,
                                                             size_t points_in,
                                                             const std::string &parameter_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), false)
{
  this->variable = variable_in;
  this->statistic = statistic_in;
  this->exact = exact_in;
  this->sweep = sweep_in;
  this->sweeps = sweeps_in;
  this->points = std::max(points_in, size_t(1));
  this->parameter = parameter_in;
}

UnnormalisedGradientEstimator::~UnnormalisedGradientEstimator()
{
}

UnnormalisedGradientEstimator::UnnormalisedGradientEstimator(const UnnormalisedGradientEstimator &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void UnnormalisedGradientEstimator::operator=(const UnnormalisedGradientEstimator &another)
{
  if (this == &another)
    return;
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* UnnormalisedGradientEstimator::duplicate() const
{
  return new UnnormalisedGradientEstimator(*this);
}

void UnnormalisedGradientEstimator::make_copy(const UnnormalisedGradientEstimator &another)
{
  this->variable = another.variable;
  this->statistic = another.statistic;
  this->exact = another.exact;
  this->sweep = another.sweep;
  this->sweeps = another.sweeps;
  this->points = another.points;
  this->parameter = another.parameter;
  this->replicates = another.replicates;
}

LikelihoodEstimatorOutput* UnnormalisedGradientEstimator::initialise()
{
  return new UnnormalisedGradientEstimatorOutput(this);
}

LikelihoodEstimatorOutput* UnnormalisedGradientEstimator::initialise(const Parameters &parameters)
{
  return this->initialise();
}

void UnnormalisedGradientEstimator::setup()
{
}

void UnnormalisedGradientEstimator::setup(const Parameters &parameters)
{
}

void UnnormalisedGradientEstimator::specific_change_data(Data* new_data)
{
}

double UnnormalisedGradientEstimator::estimate(const Parameters &parameters) const
{
  double total = 0.0;
  for (size_t m = 0; m < this->points; ++m)
  {
    arma::mat u;
    if (this->exact != NULL)
    {
      Parameters drawn = this->exact(*this->rng, parameters, *this->current_data);
      u = drawn[this->variable];
    }
    else
    {
      // sweeps (approximate): from the field in the parameters (a latent
      // field's current value), otherwise from the base measure
      Parameters state = parameters.deep_copy();
      for (size_t k = 0; k < std::max(this->sweeps, size_t(1)); ++k)
        this->sweep(*this->rng, state, *this->current_data);
      u = state[this->variable];
    }
    if (u.has_nan())
      return arma::datum::nan;
    total += this->statistic(u, parameters, *this->current_data);
  }
  return -double(this->replicates) * total / double(this->points);
}

}
