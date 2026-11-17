#include "hmm_factor_variables.h"
#include "hmm_factors.h"
#include "index.h"
#include "likelihood_estimator_output.h"
#include "particle.h"
#include "proposal_kernel.h"

namespace ilike
{
HMMFactorVariables::HMMFactorVariables()
:FactorVariables()
{
  this->likelihood_estimator_outputs.resize(0);
  this->dynamic_smcfixed_part = 0.0;
  this->hmm_factors = NULL;
}

// outputs taken by value and moved in (as VectorFactorVariables)
HMMFactorVariables::HMMFactorVariables(const HMMFactors* hmm_factors_in,
                                       std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in)
:FactorVariables()
{
  this->hmm_factors = hmm_factors_in;
  // takes ownership of the outputs
  this->likelihood_estimator_outputs.reserve(likelihood_estimator_outputs_in.size());
  for (auto i=likelihood_estimator_outputs_in.begin(); i!=likelihood_estimator_outputs_in.end(); ++i)
    this->likelihood_estimator_outputs.emplace_back(*i);
}

HMMFactorVariables::HMMFactorVariables(const HMMFactors* hmm_factors_in,
                                       std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in,
                                       Particle* particle_in)
:FactorVariables(particle_in)
{
  this->hmm_factors = hmm_factors_in;
  // takes ownership of the outputs
  this->likelihood_estimator_outputs.reserve(likelihood_estimator_outputs_in.size());
  for (auto i=likelihood_estimator_outputs_in.begin(); i!=likelihood_estimator_outputs_in.end(); ++i)
    this->likelihood_estimator_outputs.emplace_back(*i);
}

HMMFactorVariables::~HMMFactorVariables()
{
  
}

//Copy constructor for the HMMFactorVariables class.
HMMFactorVariables::HMMFactorVariables(const HMMFactorVariables &another)
:FactorVariables(another)
{
  this->make_copy(another);
}

void HMMFactorVariables::operator=(const HMMFactorVariables &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->likelihood_estimator_outputs.clear();

  FactorVariables::operator=(another);
  this->make_copy(another);
}

FactorVariables* HMMFactorVariables::copy_with_new_outputs(const Parameters &new_parameters,
                                                          const Index* new_outputs_index) const
{
  Rcpp::stop("HMMFactorVariables::copy_with_new_outputs - not written yet.");
}

FactorVariables* HMMFactorVariables::duplicate() const
{
  return( new HMMFactorVariables(*this));
}

void HMMFactorVariables::make_copy(const HMMFactorVariables &another)
{
  this->hmm_factors = another.hmm_factors;
  // a copy keeps the cached transition density, as its outputs keep their
  // cached log-likelihoods (B27: it was not copied)
  this->dynamic_smcfixed_part = another.dynamic_smcfixed_part;

  this->likelihood_estimator_outputs.resize(0);
  this->likelihood_estimator_outputs.reserve(another.likelihood_estimator_outputs.size());
  for (auto i=another.likelihood_estimator_outputs.begin();
       i!=another.likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimator_outputs.emplace_back((*i)->duplicate());
    else
      this->likelihood_estimator_outputs.emplace_back(nullptr);
  }
  
}

void HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  //if (index->size()>1)
  //  Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");

  if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
    this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                                        *this->particle->previous_self);
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(this->particle->parameters);
  }
  
}

void HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{

  if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
    this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                                  *this->particle->previous_self);
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcfixed_part(this->particle->parameters);
  }
}

double HMMFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  //if (index->size()>1)
  //  Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  double result = 0.0;

  if (this->hmm_factors->smcfixed_flag)
    result = result + this->dynamic_smcfixed_part;
  else
  {
    if ( (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
      result = result + this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                              *this->particle->previous_self);
  }
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    this->likelihood_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
    result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
  }
  
  return result;
}

double HMMFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{

  double result = 0.0;
  
  if (this->hmm_factors->smcfixed_flag)
    result = result + this->dynamic_smcfixed_part;
  else
  {
    if ( (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
      result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                        *this->particle->previous_self);
  }
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
    result = result + this->likelihood_estimator_outputs[*i]->subsample_log_likelihood;
  }
  
  return result;
}

std::vector<size_t> HMMFactorVariables::find_changed_exact_factors(const Index* new_outputs_index) const
{
  Rcpp::stop("HMMFactorVariables::find_changed_exact_factors - not written yet.");
}

double HMMFactorVariables::get_cached_likelihoods(const Index* index) const
{
  Rcpp::stop("HMMFactorVariables::get_cached_likelihoods - not written yet.");
}

double HMMFactorVariables::evaluate_likelihoods(const Index* index) const
{
  
  double result = 0.0;
  
  if ( (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
    result = result + this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                            *this->particle->previous_self);
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->likelihood_estimator_outputs[*i]->evaluate(this->particle->parameters);
  }
  
  return result;
}

double HMMFactorVariables::subsample_evaluate_likelihoods(const Index* index) const
{
  
  double result = 0.0;
  
  if ( (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
    result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                      *this->particle->previous_self);
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->likelihood_estimator_outputs[*i]->subsample_evaluate(this->particle->parameters);
  }
  
  return result;
}

arma::mat HMMFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                         const std::string &variable) const
{
  Rcpp::stop("HMMFactorVariables::direct_get_gradient_of_log - not written yet..");
  
}

arma::mat HMMFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                   const std::string &variable) const
{
  Rcpp::stop("HMMFactorVariables::direct_get_gradient_of_log - not written yet.");
  
}

bool HMMFactorVariables::recycle_all_outputs(const Factors* factors,
                                             const Parameters &new_parameters)
{
  // as HMMFactors::simulate_factor_variables, reusing the output objects
  if ((factors != this->hmm_factors) ||
      (this->likelihood_estimator_outputs.size() != this->hmm_factors->likelihood_estimators.size()))
    return false;
  for (auto i = this->likelihood_estimator_outputs.begin(); i != this->likelihood_estimator_outputs.end(); ++i)
  {
    if ((*i == NULL) || !(*i)->reinitialise(new_parameters))
      return false;
    (*i)->simulate(new_parameters);
  }
  this->dynamic_smcfixed_part = 0.0;
  return true;
}

const Factors* HMMFactorVariables::get_factors() const
{
  return this->hmm_factors;
}

void HMMFactorVariables::forget_you_were_already_written_to_file()
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    this->likelihood_estimator_outputs[i]->forget_you_were_already_written_to_file();
  }
}

void HMMFactorVariables::write_to_file(const std::string &directory_name,
                                       const std::string &index) const
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    std::string factor_directory_name = directory_name + "/factor" + toString(i+1);
    this->likelihood_estimator_outputs[i]->write_to_file(factor_directory_name,index);
  }
}

void HMMFactorVariables::close_ofstreams()
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    this->likelihood_estimator_outputs[i]->close_ofstreams();
  }
}
}
