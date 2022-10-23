#include <algorithm>
#include <numeric>
#include "smc_output.h"
#include "smc.h"
#include "importance_sampler.h"
#include "smc_mcmc_move.h"
#include "smc_marginal.h"
#include "smc_generic.h"

SMCOutput::SMCOutput()
  :LikelihoodEstimatorOutput()
{
  this->estimator = NULL;
}

SMCOutput::~SMCOutput()
{

}

SMCOutput::SMCOutput(SMC* estimator_in,
                     size_t lag_in,
                     size_t lag_proposed_in)
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_pre_last_step = 0.0;
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
  
  this->all_particles.clear();
  this->all_proposed.clear();
  //this->unnormalised_log_weights.clear();
  //this->normalised_log_weights.clear();
  //this->log_normalising_constant_ratios.clear();
  //this->incremental_log_weights.clear();

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
  //this->log_normalising_constant_ratios = another.log_normalising_constant_ratios;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  this->estimator = another.estimator;
  this->log_likelihood_pre_last_step = another.log_likelihood_pre_last_step;
}

void SMCOutput::simulate()
{
  this->estimator->simulate_smc(this);
}

void SMCOutput::evaluate_smcfixed_part()
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smc(this);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->evaluate_smcfixed_part_smc(this);
  }
}

void SMCOutput::evaluate_smcadaptive_part_given_smcfixed()
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smcadaptive_part_given_smcfixed_smc(this);
  }
  
}

void SMCOutput::simulate(const Parameters &parameters)
{
  this->estimator->simulate_smc(this, parameters);
}

void SMCOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smc(this, parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
                        //this->log_normalising_constant_ratios.end(),
                        //1.0);
  }
  else
  {
    this->estimator->evaluate_smcfixed_part_smc(this, parameters);
  }
}

void SMCOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->evaluate_smcadaptive_part_given_smcfixed_smc(this,parameters);
  }
}

void SMCOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator->subsample_simulate_smc(this, parameters);
}

void SMCOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smc(this, parameters);
    //this->log_likelihood_smcfixed_part = std::accumulate(this->log_normalising_constant_ratios.begin(),
    //this->log_normalising_constant_ratios.end(),
    //1.0);
  }
  else
  {
    this->estimator->subsample_evaluate_smcfixed_part_smc(this, parameters);
  }
}

void SMCOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(this,parameters);
  }
}

Particles* SMCOutput::add_particles()
{
  size_t num_to_pop_back = std::max<int>(0,this->all_particles.size()-lag-1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->all_particles.pop_back();
  }
  this->all_particles.push_back(Particles());
  this->all_particles.back().reserve(this->estimator->number_of_particles);
  return &this->all_particles.back();
}

void SMCOutput::add_proposed_particles(const Particles &latest_proposed_particles)
{
  size_t num_to_pop_back = std::max<int>(0,this->all_proposed.size()-lag_proposed-1);
  for (size_t i=0; i<num_to_pop_back; ++i)
  {
    this->all_proposed.pop_back();
  }
  this->all_proposed.push_back(latest_proposed_particles);
}

Particles SMCOutput::back() const
{
  return this->all_particles.back();
}

Particles& SMCOutput::back()
{
  return this->all_particles.back();
}

std::deque<Particles>::iterator SMCOutput::end()
{
  return this->all_particles.end();
}

std::deque<Particles>::const_iterator SMCOutput::end() const
{
  return this->all_particles.end();
}

double SMCOutput::latest_log_normalising_constant_ratio() const
{
  return this->all_particles.back().log_normalising_constant_ratio;
}

void SMCOutput::update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
{
  this->all_particles.back().update_weights(latest_unnormalised_log_incremental_weights);
}

//void SMCOutput::initialise_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
//{
  //arma::colvec latest_unnormalised_log_weights;
  //if (this->unnormalised_log_incremental_weights.size()>0)
  //  latest_unnormalised_log_weights = this->unnormalised_log_weights.back() + latest_unnormalised_log_incremental_weights;
  //else
  //  latest_unnormalised_log_weights = latest_unnormalised_log_weight_updates;
  //this->unnormalised_log_incremental_weights.push_back(latest_unnormalised_log_incremental_weights);
  //size_t num_to_pop_back = std::max<int>(0,unnormalised_log_incremental_weights.size()-lag);
  //for (size_t i=0; i<num_to_pop_back; ++i)
  //{
  //  this->unnormalised_log_incremental_weights.pop_back();
  //}
//}

//void SMCOutput::set_unnormalised_log_incremental_weights(const arma::colvec &latest_unnormalised_log_incremental_weights)
//{
//  this->unnormalised_log_incremental_weights.push_back(latest_unnormalised_log_incremental_weights);
//  size_t num_to_pop_back = std::max<int>(0,this->unnormalised_log_incremental_weights.size()-lag);
//  for (size_t i=0; i<num_to_pop_back; ++i)
//  {
//    this->unnormalised_log_incremental_weights.pop_back();
//  }
//}

//void SMCOutput::initialise_next_step()
//{
//  arma::colvec init(this->estimator->number_of_particles);
//  init.fill(0.0);
//  this->unnormalised_log_weights.push_back(init);
//  this->incremental_log_weights.push_back(init);
//  //this->unnormalised_log_incremental_weights.push_back(init);
//  this->log_normalising_constant_ratios.push_back(0.0);
//}

void SMCOutput::normalise_weights()
{
  this->log_likelihood_pre_last_step = this->log_likelihood;
  this->all_particles.back().normalise_weights();
}

void SMCOutput::resample()
{
  this->estimator->resample(this);
}

LikelihoodEstimator* SMCOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

size_t SMCOutput::number_of_smc_iterations() const
{
  return this->all_particles.size();
}

arma::mat SMCOutput::get_gradient_of_log(const std::string &variable,
                                         const Parameters &x)
{
  throw std::runtime_error("SMCOutput::get_gradient_of_log - not yet implemented.");
}

arma::mat SMCOutput::subsample_get_gradient_of_log(const std::string &variable,
                                         const Parameters &x)
{
  throw std::runtime_error("SMCOutput::get_gradient_of_log - not yet implemented.");
}

void SMCOutput::print(std::ostream &os) const
{
  os << "all_particles" << std::endl << "(" << std::endl;
  std::deque<Particles>::const_iterator it;
  for (it=this->all_particles.begin();it!=this->all_particles.end();++it)
  {
    if (it==this->all_particles.begin())
      os << *it;
    else
      os << std::endl << "," << std::endl << *it;
  }
  os << std::endl << ")" << std::endl;

  os << "all_proposed" << std::endl << "(" << std::endl;
  for (it=this->all_proposed.begin();it!=this->all_proposed.end();++it)
  {
    if (it==this->all_proposed.begin())
      os << *it;
    else
      os << std::endl << "," << std::endl << *it;
  }
  os << std::endl << ")" << std::endl;

  /*
  os << "unnormalised_log_weights" << std::endl << "(" << std::endl;
  std::deque<arma::colvec>::const_iterator itd;
  for (itd=this->unnormalised_log_weights.begin();itd!=this->unnormalised_log_weights.end();++itd)
  {
    if (itd==this->unnormalised_log_weights.begin())
      os << *itd;
    else
      os << std::endl << "," << std::endl << *itd;
  }
  os << std::endl << ")" << std::endl;

  os << "normalised_log_weights" << std::endl << "(" << std::endl;
  for (itd=this->normalised_log_weights.begin();itd!=this->normalised_log_weights.end();++itd)
  {
    if (itd==this->normalised_log_weights.begin())
      os << *itd;
    else
      os << std::endl << "," << std::endl << *itd;
  }
  os << std::endl << ")" << std::endl;

  os << "log_normalising_constant_ratios" << std::endl << "(" << std::endl;
  std::vector<double>::const_iterator i;
  for (i=this->log_normalising_constant_ratios.begin();i!=this->log_normalising_constant_ratios.end();++i)
  {
    if (i==this->log_normalising_constant_ratios.begin())
      os << *i;
    else
      os << std::endl << "," << std::endl << *i;
  }
  os << std::endl << ")" << std::endl;
  */
}
