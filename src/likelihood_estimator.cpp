#include "likelihood_estimator.h"
//#include "model_and_algorithm.h"
#include "likelihood_estimator_output.h"
#include "smc_worker.h"
#include "factors.h"

LikelihoodEstimator::LikelihoodEstimator()
{
  this->factors = NULL;
}

LikelihoodEstimator::LikelihoodEstimator(RandomNumberGenerator* rng_in,
                                         size_t* seed_in,
                                         Data* data_in,
                                         bool smcfixed_flag_in)
{
  this->data = data_in;
  this->current_data = this->data;
  //this->model_and_algorithm = ModelAndAlgorithm() ;
  this->rng = rng_in;
  this->seed = seed_in;
  this->smcfixed_flag = smcfixed_flag_in;
  this->subsampler = NULL;
  this->factors = NULL;
}

LikelihoodEstimator::~LikelihoodEstimator()
{
  if (this->factors!=NULL)
    delete this->factors;
}

LikelihoodEstimator::LikelihoodEstimator(const LikelihoodEstimator &another)
{
  this->make_copy(another);
}

void LikelihoodEstimator::operator=(const LikelihoodEstimator &another)
{
  if(this == &another)
    return;
  
  if (this->factors!=NULL)
    delete this->factors;

  this->make_copy(another);
}

void LikelihoodEstimator::make_copy(const LikelihoodEstimator &another)
{
  this->data = another.data;
  this->current_data = another.current_data;
  //this->model_and_algorithm = another.model_and_algorithm;
  this->rng = another.rng;
  this->seed = another.seed;
  this->subsampler = another.subsampler;
  this->smcfixed_flag = another.smcfixed_flag;
  if (another.factors!=NULL)
    this->factors = another.factors->duplicate();
  else
    this->factors = NULL;
}

/*
double LikelihoodEstimator::estimate()
{
  LikelihoodEstimatorOutput* output = this->initialise();
  output->simulate();
  output->evaluate();
  double log_llhd = output->log_likelihood;
  delete output;
  return log_llhd;
}
 */

double LikelihoodEstimator::estimate(const Parameters &parameters)
{
  LikelihoodEstimatorOutput* output = this->initialise(parameters);
  output->simulate(parameters);
  output->evaluate(parameters);
  double log_llhd = output->log_likelihood;
  delete output;
  return log_llhd;
}

void LikelihoodEstimator::change_data()
{
  this->current_data = this->data;
}

void LikelihoodEstimator::change_data(Data* new_data)
{
  this->current_data = new_data;
}

Data* LikelihoodEstimator::get_data() const
{
  return this->data;
}

bool LikelihoodEstimator::get_smcfixed_flag() const
{
  return this->smcfixed_flag;
}
