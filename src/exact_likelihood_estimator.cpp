#include "exact_likelihood_estimator.h"
#include "exact_likelihood_estimator_output.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator()
  :LikelihoodEstimator()
{
}

ExactLikelihoodEstimator::ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                                                   size_t* seed_in,
                                                   const Data* data_in,
                                                   EvaluateLogLikelihoodPtr func_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->func = func_in;
  //this->output = new ExactLikelihoodEstimatorOutput();
}

ExactLikelihoodEstimator::~ExactLikelihoodEstimator()
{
}

//Copy constructor for the ExactLikelihoodEstimator class.
ExactLikelihoodEstimator::ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void ExactLikelihoodEstimator::operator=(const ExactLikelihoodEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* ExactLikelihoodEstimator::duplicate(void)const
{
  return( new ExactLikelihoodEstimator(*this));
}

void ExactLikelihoodEstimator::make_copy(const ExactLikelihoodEstimator &another)
{
  this->func = another.func;
  //if (this->output!=NULL)
  //  this->output = another.output->duplicate();
}

// double ExactLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* ExactLikelihoodEstimator::initial_simulate(const Parameters &parameters)
{
  return new ExactLikelihoodEstimatorOutput(this);
}

// void ExactLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
