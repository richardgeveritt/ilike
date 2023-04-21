#include "annealed_likelihood_estimator.h"
#include "annealed_likelihood_estimator_output.h"

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
:LikelihoodEstimator(rng_in, seed_in, data_in, smcfixed_flag_in)
{
  this->estimator = estimator_in;
  this->function_power = function_power_in;
  this->power_variable = power_variable_in;
  this->constant_power = 1.0;
  this->use_constant = FALSE;
  //this->smcfixed_flag = this->estimator->get_smcfixed_flag();
  //this->output = new AnnealedLikelihoodEstimatorOutput();
}

AnnealedLikelihoodEstimator::AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                         size_t* seed_in,
                                                         Data* data_in,
                                                         LikelihoodEstimator* estimator_in,
                                                         double constant_power_in,
                                                         bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, smcfixed_flag_in)
{
  this->estimator = estimator_in;
  this->function_power = PowerFunctionPtr();
  this->power_variable = "";
  this->constant_power = constant_power_in;
  this->use_constant = FALSE;
  //this->smcfixed_flag = this->estimator->get_smcfixed_flag();
  //this->output = new AnnealedLikelihoodEstimatorOutput();
}

AnnealedLikelihoodEstimator::~AnnealedLikelihoodEstimator()
{
  if (this->estimator!=NULL)
    delete this->estimator;
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
  
  if (this->estimator!=NULL)
    delete this->estimator;

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
    this->estimator = another.estimator->duplicate();
  else
    this->estimator = NULL;
  
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

// void AnnealedLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
