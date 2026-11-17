#include "annealed_likelihood_estimator.h"
#include "annealed_likelihood_estimator_output.h"

namespace ilike
{
AnnealedLikelihoodEstimator::AnnealedLikelihoodEstimator()
:LikelihoodEstimator()
{
}

AnnealedLikelihoodEstimator::AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                         size_t* seed_in,
                                                         Data* data_in,
                                                         LikelihoodEstimator* estimator_in,
                                                         PowerFunctionPtr function_power_in,
                                                         const std::string &power_variable_in,
                                                         bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->estimator.reset(estimator_in);
  this->function_power = function_power_in;
  this->power_variable = power_variable_in;
  this->constant_power = 1.0;
  this->use_constant = FALSE;
}

AnnealedLikelihoodEstimator::AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                         size_t* seed_in,
                                                         Data* data_in,
                                                         LikelihoodEstimator* estimator_in,
                                                         double constant_power_in,
                                                         bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->estimator.reset(estimator_in);
  this->function_power = PowerFunctionPtr();
  this->power_variable = "";
  this->constant_power = constant_power_in;
  this->use_constant = FALSE;
}

AnnealedLikelihoodEstimator::~AnnealedLikelihoodEstimator()
{
  this->estimator.reset();
}

//Copy constructor for the AnnealedLikelihoodEstimator class.
AnnealedLikelihoodEstimator::AnnealedLikelihoodEstimator(const AnnealedLikelihoodEstimator &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void AnnealedLikelihoodEstimator::operator=(const AnnealedLikelihoodEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->estimator.reset();
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* AnnealedLikelihoodEstimator::duplicate() const
{
  return( new AnnealedLikelihoodEstimator(*this));
}

void AnnealedLikelihoodEstimator::make_copy(const AnnealedLikelihoodEstimator &another)
{
  if (another.estimator!=NULL)
    this->estimator.reset(another.estimator->duplicate());
  else
    this->estimator.reset();
  
  this->function_power = another.function_power;
  this->power_variable = another.power_variable;
  this->constant_power = another.constant_power;
  this->use_constant = another.use_constant;
}

// double AnnealedLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* AnnealedLikelihoodEstimator::initialise()
{
  LikelihoodEstimatorOutput* estimator_output = this->estimator->initialise();
  return new AnnealedLikelihoodEstimatorOutput(this,estimator_output);
}

LikelihoodEstimatorOutput* AnnealedLikelihoodEstimator::initialise(const Parameters &parameters)
{
  LikelihoodEstimatorOutput* estimator_output = this->estimator->initialise(parameters);
  return new AnnealedLikelihoodEstimatorOutput(this,estimator_output);
}

void AnnealedLikelihoodEstimator::setup()
{
  this->estimator->setup();
}

void AnnealedLikelihoodEstimator::setup(const Parameters &parameters)
{
  this->estimator->setup(parameters);
}

void AnnealedLikelihoodEstimator::specific_change_data(Data* new_data)
{
  this->estimator->change_data_with_raw_pointer(new_data);
}

bool AnnealedLikelihoodEstimator::supports_restore_from_summary() const
{
  // the summary is the inner output's value (inner smcfixed) or the inner
  // output's own summary; either way restoring the inner output restores it
  return (this->estimator!=NULL) &&
    this->estimator->supports_restore_from_summary();
}

bool AnnealedLikelihoodEstimator::supports_parallel_chunks() const
{
  return (this->estimator!=NULL) && this->estimator->supports_parallel_chunks();
}

bool AnnealedLikelihoodEstimator::supports_parallel_copies() const
{
  return (this->estimator!=NULL) && this->estimator->supports_parallel_copies();
}

void AnnealedLikelihoodEstimator::redirect_random_numbers(RandomNumberGenerator* rng_in,
                                                          size_t* seed_in)
{
  LikelihoodEstimator::redirect_random_numbers(rng_in, seed_in);
  if (this->estimator!=NULL)
    this->estimator->redirect_random_numbers(rng_in, seed_in);
}

bool AnnealedLikelihoodEstimator::supports_flat_generations() const
{
  // its smcfixed part is the inner estimator's (evaluated in full when the
  // inner one is smcfixed), and its summary the inner one's
  return (this->estimator!=NULL) && this->estimator->supports_flat_generations();
}

}

