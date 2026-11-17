#include "hmm_factors.h"
#include "proposal_kernel.h"
#include "hmm_factor_variables.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

namespace ilike
{
HMMFactors::HMMFactors()
:Factors()
{
  this->likelihood_estimators.resize(0);
  this->likelihood_estimator_temp_data.resize(0);
}

HMMFactors::HMMFactors(ProposalKernel* transition_kernel_in,
                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:Factors()
{
  this->transition_kernel.reset(transition_kernel_in);
  // takes ownership of the estimators
  this->likelihood_estimators.reserve(likelihood_estimators_in.size());
  for (auto i=likelihood_estimators_in.begin(); i!=likelihood_estimators_in.end(); ++i)
    this->likelihood_estimators.emplace_back(*i);
  this->likelihood_estimator_temp_data.resize(0);
}

HMMFactors::~HMMFactors()
{
  this->transition_kernel.reset();
  
}

//Copy constructor for the HMMFactors class.
HMMFactors::HMMFactors(const HMMFactors &another)
:Factors(another)
{
  this->make_copy(another);
}

void HMMFactors::operator=(const HMMFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->transition_kernel.reset();
  this->likelihood_estimators.clear();

  Factors::operator=(another);
  this->make_copy(another);
}

Factors* HMMFactors::duplicate() const
{
  return( new HMMFactors(*this));
}

void HMMFactors::make_copy(const HMMFactors &another)
{
  if (another.transition_kernel!=NULL)
    this->transition_kernel.reset(another.transition_kernel->proposal_kernel_duplicate());
  else
    this->transition_kernel.reset();
  
  this->smcfixed_flag = another.smcfixed_flag;

  this->likelihood_estimators.resize(0);
  this->likelihood_estimators.reserve(another.likelihood_estimators.size());
  for (auto i=another.likelihood_estimators.begin();
       i!=another.likelihood_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimators.emplace_back((*i)->duplicate());
    else
      this->likelihood_estimators.emplace_back(nullptr);
  }
  
  this->likelihood_estimator_temp_data = another.likelihood_estimator_temp_data;
}

void HMMFactors::set_data(const Index* index)
{
  if (this->likelihood_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  
  // assumption that all llhd_estimators point to the same data
  Data* all_data = this->likelihood_estimators[0]->get_data();
  
  if (this->likelihood_estimator_temp_data.size()==0)
  {
    this->likelihood_estimator_temp_data.push_back(std::make_shared<Data>(all_data->rows(indices)));
  }
  else
  {
    this->likelihood_estimator_temp_data[0] = std::make_shared<Data>(all_data->rows(indices));
  }
  
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->change_data(this->likelihood_estimator_temp_data[0]);
  }
}

FactorVariables* HMMFactors::simulate_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  for (auto i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->simulate(simulated_parameters);
    outputs.back()->write_to_file_flag = false;
  }

  return new HMMFactorVariables(this,
                                std::move(outputs));
}

FactorVariables* HMMFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  
  for (auto i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->subsample_simulate(simulated_parameters);
    outputs.back()->write_to_file_flag = false;
  }

  return new HMMFactorVariables(this,
                                std::move(outputs));
}

void HMMFactors::setup()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void HMMFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}

void HMMFactors::specific_change_data(Data* new_data)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->change_data_with_raw_pointer(new_data);
  }
}

Data* HMMFactors::get_current_data()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    return (*i)->get_current_data();
  }
  
  return NULL;
}

bool HMMFactors::supports_parallel_chunks() const
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if ((*i==NULL) || !(*i)->supports_parallel_chunks())
      return false;
  }
  return true;
}

bool HMMFactors::supports_parallel_copies() const
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if ((*i==NULL) || !(*i)->supports_parallel_copies())
      return false;
  }
  return true;
}

void HMMFactors::redirect_random_numbers(RandomNumberGenerator* rng_in,
                                         size_t* seed_in)
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if (*i!=NULL)
      (*i)->redirect_random_numbers(rng_in, seed_in);
  }
}

}
