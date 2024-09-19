#include "abc_likelihood_estimator.h"
#include "abc_likelihood_estimator_output.h"

namespace ilike
{

ABCLikelihoodEstimator::ABCLikelihoodEstimator()
:LikelihoodEstimator()
{
}

ABCLikelihoodEstimator::ABCLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               ABCKernelFactor* abc_kernel_in,
                                               bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->abc_kernel = abc_kernel_in;
}

ABCLikelihoodEstimator::~ABCLikelihoodEstimator()
{
  if (this->abc_kernel!=NULL)
    delete this->abc_kernel;
}

//Copy constructor for the ABCLikelihoodEstimator class.
ABCLikelihoodEstimator::ABCLikelihoodEstimator(const ABCLikelihoodEstimator &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void ABCLikelihoodEstimator::operator=(const ABCLikelihoodEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->abc_kernel!=NULL)
    delete this->abc_kernel;
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* ABCLikelihoodEstimator::duplicate() const
{
  return( new ABCLikelihoodEstimator(*this));
}

void ABCLikelihoodEstimator::make_copy(const ABCLikelihoodEstimator &another)
{
  if (another.abc_kernel!=NULL)
    this->abc_kernel = another.abc_kernel->abc_kernel_factor_duplicate();
  else
    this->abc_kernel = NULL;
}

// double ABCLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* ABCLikelihoodEstimator::initialise()
{
  return new ABCLikelihoodEstimatorOutput(this);
}

LikelihoodEstimatorOutput* ABCLikelihoodEstimator::initialise(const Parameters &parameters)
{
  return new ABCLikelihoodEstimatorOutput(this);
}

void ABCLikelihoodEstimator::setup()
{
  
}

void ABCLikelihoodEstimator::setup(const Parameters &parameters)
{
  
}

void ABCLikelihoodEstimator::specific_change_data(Data* new_data)
{
  this->abc_kernel->set_data(new_data);
}

}
