#include <algorithm>
#include "smc_output.h"
#include "utils.h"
#include "smc.h"

SMCOutput::SMCOutput()
  :LikelihoodEstimatorOutput()
{
}

SMCOutput::~SMCOutput()
{

}

SMCOutput::SMCOutput(SMC* estimator_in,
                     size_t lag_in,
                     size_t lag_proposed_in)
  :LikelihoodEstimatorOutput()
{
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->estimator = estimator_in;
}

//Copy constructor for the SMCOutput class.
SMCOutput::SMCOutput(const SMCOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void SMCOutput::operator=(const SMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* SMCOutput::duplicate(void)const
{
  return( new SMCOutput(*this));
}

SMCOutput* SMCOutput::smc_duplicate(void)const
{
  return( new SMCOutput(*this));
}

void SMCOutput::make_copy(const SMCOutput &another)
{
  this->all_particles = another.all_particles;
  this->all_proposed = another.all_proposed;
  this->unnormalised_log_weights = another.unnormalised_log_weights;
  this->normalised_log_weights = another.normalised_log_weights;
  this->log_normalising_constant_ratios = another.log_normalising_constant_ratios;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  this->estimator = another.estimator;
}

void SMCOutput::continue_simulate(const Parameters &parameters)
{

}

void SMCOutput::estimate(const Parameters &parameters)
{
}

void SMCOutput::add_particles(const Particles &latest_particles)
{
  this->all_particles.push_back(latest_particles);
  size_t num_to_pop_back = std::max<int>(0,all_particles.size()-lag);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->all_particles.pop_back();
  }
}

void SMCOutput::add_proposed_particles(const Particles &latest_proposed_particles)
{
  this->all_proposed.push_back(latest_proposed_particles);
  size_t num_to_pop_back = std::max<int>(0,all_proposed.size()-lag_proposed);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->all_proposed.pop_back();
  }
}

void SMCOutput::add_weights(const arma::colvec &latest_unnormalised_log_weight_updates)
{
  arma::colvec latest_unnormalised_log_weights;
  if (this->unnormalised_log_weights.size()>0)
    latest_unnormalised_log_weights = this->unnormalised_log_weights.back() + latest_unnormalised_log_weight_updates;
  else
    latest_unnormalised_log_weights = latest_unnormalised_log_weight_updates;
  this->unnormalised_log_weights.push_back(latest_unnormalised_log_weights);
  size_t num_to_pop_back = std::max<int>(0,unnormalised_log_weights.size()-lag);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->unnormalised_log_weights.pop_back();
  }

  double log_normalising_constant_ratio = log_sum_exp(latest_unnormalised_log_weights);
  this->log_normalising_constant_ratios.push_back(log_normalising_constant_ratio);
  arma::colvec latest_normalised_log_weights = this->unnormalised_log_weights.back() - log_normalising_constant_ratio;
  this->normalised_log_weights.push_back(latest_normalised_log_weights);
  num_to_pop_back = std::max<int>(0,normalised_log_weights.size()-lag);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->normalised_log_weights.pop_back();
  }
}
