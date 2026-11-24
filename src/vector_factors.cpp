#include "vector_factors.h"
#include "likelihood_estimator.h"
#include "vector_factor_variables.h"
#include "index.h"
#include "likelihood_estimator_output.h"

namespace ilike
{
VectorFactors::VectorFactors()
:Factors()
{
  this->likelihood_estimators.resize(0);
  this->likelihood_estimator_temp_data.resize(0);
}

VectorFactors::VectorFactors(const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:Factors()
{
  // takes ownership of the estimators
  this->likelihood_estimators.reserve(likelihood_estimators_in.size());
  for (auto i=likelihood_estimators_in.begin(); i!=likelihood_estimators_in.end(); ++i)
    this->likelihood_estimators.emplace_back(*i);
  this->find_deferring();
}

void VectorFactors::find_deferring()
{
  this->any_defers = false;
  for (const auto &estimator : this->likelihood_estimators)
    if ((estimator != NULL) && estimator->defers_simulation())
      this->any_defers = true;
}

VectorFactors::~VectorFactors()
{
  
}

//Copy constructor for the VectorFactors class.
VectorFactors::VectorFactors(const VectorFactors &another)
:Factors(another)
{
  this->make_copy(another);
}

void VectorFactors::operator=(const VectorFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->likelihood_estimators.clear();

  Factors::operator=(another);
  this->make_copy(another);
}

Factors* VectorFactors::duplicate() const
{
  return( new VectorFactors(*this));
}

void VectorFactors::make_copy(const VectorFactors &another)
{
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
  this->find_deferring();
  
  this->likelihood_estimator_temp_data = another.likelihood_estimator_temp_data;
}

void VectorFactors::set_data(const Index* index)
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

LikelihoodEstimator* VectorFactors::get_likelihood_estimator(size_t i) const
{
  return this->likelihood_estimators[i].get();
}

bool VectorFactors::supports_parallel_copies() const
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if ((*i==NULL) || !(*i)->supports_parallel_copies())
      return false;
  }
  return true;
}

void VectorFactors::redirect_random_numbers(RandomNumberGenerator* rng_in,
                                            size_t* seed_in)
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if (*i!=NULL)
      (*i)->redirect_random_numbers(rng_in, seed_in);
  }
}

bool VectorFactors::supports_parallel_chunks() const
{
  for (auto i=this->likelihood_estimators.begin(); i!=this->likelihood_estimators.end(); ++i)
  {
    if ((*i==NULL) || !(*i)->supports_parallel_chunks())
      return false;
  }
  return true;
}

bool VectorFactors::supports_flat_moves(const Index* index, const Index* other_index) const
{
  if (index==NULL)
    return false;
  std::vector<int> role(this->likelihood_estimators.size(), 0);
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    if ((*i >= role.size()) || (this->likelihood_estimators[*i]==NULL) ||
        !this->likelihood_estimators[*i]->supports_restore_from_summary())
      return false;
    role[*i] = 1;
  }
  if (other_index!=NULL)
  {
    for (auto i=other_index->begin(); i!=other_index->end(); ++i)
    {
      if ((*i >= role.size()) || (role[*i]!=0) || (this->likelihood_estimators[*i]==NULL) ||
          !this->likelihood_estimators[*i]->supports_reinitialise())
        return false;
      role[*i] = 2;
    }
  }
  for (auto r=role.begin(); r!=role.end(); ++r)
  {
    if (*r==0)
      return false;
  }
  return true;
}

bool VectorFactors::supports_flat_generations(const Index* index) const
{
  if (index==NULL)
    return false;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    if ((*i >= this->likelihood_estimators.size()) || (this->likelihood_estimators[*i]==NULL) ||
        !this->likelihood_estimators[*i]->supports_flat_generations())
      return false;
  }
  return true;
}

FactorVariables* VectorFactors::simulate_factor_variables(const Parameters &simulated_parameters) const
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

  return new VectorFactorVariables(this,
                                   std::move(outputs));
}

FactorVariables* VectorFactors::subsample_simulate_factor_variables(const Parameters &simulated_parameters) const
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

  return new VectorFactorVariables(this,
                                   std::move(outputs));
}

void VectorFactors::setup()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void VectorFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}

void VectorFactors::specific_change_data(Data* new_data)
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    (*i)->change_data_with_raw_pointer(new_data);
  }
}

Data* VectorFactors::get_current_data()
{
  for (auto i=this->likelihood_estimators.begin();
       i!=this->likelihood_estimators.end();
       ++i)
  {
    return (*i)->get_current_data();
  }
  
  return NULL;
}

double VectorFactors::log_maximum(const Index* index, const Parameters &parameters) const
{
  double output = 0.0;
  for (auto i = index->begin(); i != index->end(); ++i)
  {
    if ((*i >= this->likelihood_estimators.size()) || (this->likelihood_estimators[*i] == NULL))
      return arma::datum::inf;
    output = output + this->likelihood_estimators[*i]->log_maximum(parameters);
  }
  return output;
}

// (dev/plan.md 3.21 (c2)) a factor's estimator keeps an inner state following
// the SMC target (checked once per SMC step, not per particle)
bool VectorFactors::advances_at_caller_target() const
{
  for (const auto &estimator : this->likelihood_estimators)
    if (estimator->advances_at_caller_target())
      return true;
  return false;
}


bool VectorFactors::any_estimate() const
{
  for (const auto &estimator : this->likelihood_estimators)
    if (estimator->is_estimate())
      return true;
  return false;
}


void VectorFactors::set_outputs_written(bool written)
{
  for (auto &estimator : this->likelihood_estimators)
    if (estimator != NULL)
      estimator->set_outputs_written(written);
}

}
