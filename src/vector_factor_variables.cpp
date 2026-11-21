#include <algorithm>
#include "vector_factor_variables.h"
#include "likelihood_estimator_output.h"
#include "index.h"
#include "particle.h"
#include "vector_factors.h"
#include "likelihood_estimator.h"
#include "exact_likelihood_estimator_output.h"

namespace ilike
{
VectorFactorVariables::VectorFactorVariables()
:FactorVariables()
{
  this->likelihood_estimator_outputs.resize(0);
  this->vector_factors = NULL;
}

VectorFactorVariables::~VectorFactorVariables()
{
}

// The outputs vector is taken by value and moved in: callers that pass a
// temporary or std::move their vector no longer have it copied.
VectorFactorVariables::VectorFactorVariables(const VectorFactors* vector_factors_in,
                                             std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in)
:FactorVariables()
{
  this->vector_factors = vector_factors_in;
  // takes ownership of the outputs
  this->likelihood_estimator_outputs.reserve(likelihood_estimator_outputs_in.size());
  for (auto i=likelihood_estimator_outputs_in.begin(); i!=likelihood_estimator_outputs_in.end(); ++i)
    this->likelihood_estimator_outputs.emplace_back(*i);
}

VectorFactorVariables::VectorFactorVariables(VectorFactors* vector_factors_in,
                                             std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in,
                                             Particle* particle_in)
:FactorVariables(particle_in)
{
  this->vector_factors = vector_factors_in;
  // takes ownership of the outputs
  this->likelihood_estimator_outputs.reserve(likelihood_estimator_outputs_in.size());
  for (auto i=likelihood_estimator_outputs_in.begin(); i!=likelihood_estimator_outputs_in.end(); ++i)
    this->likelihood_estimator_outputs.emplace_back(*i);
}

//Copy constructor for the VectorFactorVariables class.
VectorFactorVariables::VectorFactorVariables(const VectorFactorVariables &another)
:FactorVariables(another)
{
  this->make_copy(another);
}

void VectorFactorVariables::operator=(const VectorFactorVariables &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->likelihood_estimator_outputs.clear();
  
  FactorVariables::operator=(another);
  this->make_copy(another);
}

// The variables whose matrices differ between the current point and a
// proposed one (C205). A proposed point is the current one with the proposed
// variables written over: the others share the current point's matrices, so
// comparing pointers finds the changed ones without looking any name up.
// False when the two do not have the same variables in the same order (then
// each factor compares every variable it reads, as before).
static bool changed_variables(const Parameters &before,
                              const Parameters &after,
                              std::vector<std::string> &changed)
{
  changed.clear();
  if (before.vector_size() != after.vector_size())
    return false;
  auto b = before.vector_begin();
  for (auto a = after.vector_begin(); a != after.vector_end(); ++a, ++b)
  {
    if (a->first != b->first)
      return false;
    if (a->second.first.get() != b->second.first.get())
      changed.push_back(a->first);
  }
  return true;
}

FactorVariables* VectorFactorVariables::copy_with_new_outputs(const Parameters &new_parameters,
                                                             const Index* new_outputs_index) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimator_outputs.size());
  std::vector<std::string> changed;
  int changed_known = -1; // (found when first needed)
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    bool new_output = std::find(new_outputs_index->begin(),
                                new_outputs_index->end(),
                                i)!=new_outputs_index->end();
    if (new_output)
    {
      outputs.push_back(this->vector_factors->get_likelihood_estimator(i)->initialise(new_parameters));
      outputs.back()->simulate(new_parameters);
      outputs.back()->write_to_file_flag = false;
      if ((this->likelihood_estimator_outputs[i]!=NULL) && (this->particle!=NULL) && outputs.back()->carries_over())
      {
        if (changed_known < 0)
          changed_known = changed_variables(this->particle->parameters, new_parameters, changed) ? 1 : 0;
        outputs.back()->carry_over_from(*this->likelihood_estimator_outputs[i],
                                        this->particle->parameters,
                                        new_parameters,
                                        (changed_known == 1) ? &changed : NULL);
      }
    }
    else if (this->likelihood_estimator_outputs[i]!=NULL)
    {
      outputs.push_back(this->likelihood_estimator_outputs[i]->duplicate());
    }
    else
    {
      outputs.push_back(NULL);
    }
  }
  return new VectorFactorVariables(this->vector_factors,
                                   std::move(outputs));
}

FactorVariables* VectorFactorVariables::duplicate() const
{
  return( new VectorFactorVariables(*this));
}

void VectorFactorVariables::make_copy(const VectorFactorVariables &another)
{
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
  this->vector_factors = another.vector_factors;
}

void VectorFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (index==NULL)
  {
  }
  else
  {
    for (auto i=index->begin();
         i!=index->end();
         ++i)
    {
      if (this->likelihood_estimator_outputs[*i]!=NULL)
      {
        this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(this->particle->parameters);
      }
    }
  }
}

void VectorFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcfixed_part(this->particle->parameters);
    }
  }
}

double VectorFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    
    this->likelihood_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
    result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
    
  }
  return result;
}

double VectorFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
      result = result + this->likelihood_estimator_outputs[*i]->subsample_log_likelihood;
    }
  }
  return result;
}

std::vector<size_t> VectorFactorVariables::find_changed_exact_factors(const Index* new_outputs_index) const
{
  std::vector<size_t> changed;
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    if (std::find(new_outputs_index->begin(),
                  new_outputs_index->end(),
                  i)!=new_outputs_index->end())
      continue;
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[i].get();
    if (dynamic_cast<ExactLikelihoodEstimatorOutput*>(output)==NULL)
      continue;
    LikelihoodEstimatorOutput* fresh = output->duplicate();
    double fresh_value = fresh->evaluate(this->particle->parameters);
    delete fresh;
    double cached_value = output->log_likelihood;
    if ( (fresh_value!=cached_value) &&
         !(std::abs(fresh_value-cached_value) <= 1e-8*(1.0+std::abs(fresh_value))) )
      changed.push_back(i);
  }
  return changed;
}

double VectorFactorVariables::get_cached_likelihoods(const Index* index) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
    }
  }
  return result;
}

double VectorFactorVariables::evaluate_likelihoods(const Index* index) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      result = result + this->likelihood_estimator_outputs[*i]->evaluate(this->particle->parameters);
    }
  }
  return result;
}

bool VectorFactorVariables::simulation_pending(const Index* index) const
{
  if ((this->vector_factors == NULL) || !this->vector_factors->can_defer())
    return false; // (no estimator defers its draws: one flag, not a loop)
  for (auto i=index->begin(); i!=index->end(); ++i)
    if ((this->likelihood_estimator_outputs[*i]!=NULL) && this->likelihood_estimator_outputs[*i]->simulation_pending())
      return true;
  return false;
}

// (dev/plan.md 3.21) early rejection: the outputs that have made their draws
// (priors, exact factors, plates) first, in index order; then, while the sum
// plus the upper bounds of the outputs not yet evaluated can still reach the
// threshold, the outputs with deferred draws, each told how low it may go
double VectorFactorVariables::evaluate_likelihoods_with_threshold(const Index* index, double threshold, bool &below) const
{
  below = false;
  double result = 0.0;
  std::vector<LikelihoodEstimatorOutput*> pending;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    if (output == NULL)
      continue;
    if (output->simulation_pending())
      pending.push_back(output);
    else
      result = result + output->evaluate(this->particle->parameters);
  }
  if (pending.empty())
    return result; // (the same sum, in the same order, as evaluate_likelihoods)
  if (result == -arma::datum::inf)
  {
    below = true; // (e.g. a prior that is zero: nothing more is needed)
    return -arma::datum::inf;
  }
  // later[k]: the sum of the upper bounds of the pending outputs after k
  // (summed, not found by subtraction: inf - inf would be meaningless)
  std::vector<double> later(pending.size() + 1, 0.0);
  for (size_t k = pending.size(); k > 0; --k)
    later[k - 1] = later[k] + pending[k - 1]->log_upper_bound(this->particle->parameters);
  for (size_t k = 0; k < pending.size(); ++k)
  {
    if (result + later[k] < threshold)
    {
      below = true;
      return -arma::datum::inf;
    }
    double remaining = later[k + 1];
    bool output_below = false;
    double value = pending[k]->evaluate_with_threshold(this->particle->parameters, threshold - result - remaining, output_below);
    if (output_below)
    {
      below = true;
      return -arma::datum::inf;
    }
    result = result + value;
  }
  if (result < threshold)
    below = true;
  return result;
}

double VectorFactorVariables::subsample_evaluate_likelihoods(const Index* index) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      result = result + this->likelihood_estimator_outputs[*i]->subsample_evaluate(this->particle->parameters);
    }
  }
  return result;
}

arma::mat VectorFactorVariables::direct_get_gradient_of_log(const std::string &variable) const
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->get_gradient_of_log(variable,
                                                  this->particle->parameters);
    }
  }
  return result;
}

arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable) const
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
    {
      result = result + (*i)->subsample_get_gradient_of_log(variable,
                                                            this->particle->parameters);
    }
  }
  return result;
}

arma::mat VectorFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                            const std::string &variable) const
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      result = result + this->likelihood_estimator_outputs[*i]->get_gradient_of_log(variable,
                                                                                    this->particle->parameters);
    }
  }
  return result;
}

// several variables at once (C204): each factor in the index is asked once
// for all of them, and the gradients are summed over the factors
std::vector<arma::mat> VectorFactorVariables::direct_get_gradients_of_log(const Index* index,
                                                                         const std::vector<std::string> &variables) const
{
  return this->direct_get_gradients_of_log_at(index, variables, this->particle->parameters);
}

// at a given point (the factors' gradients are functions of the parameters
// they are given; C208)
std::vector<arma::mat> VectorFactorVariables::direct_get_gradients_of_log_at(const Index* index,
                                                                            const std::vector<std::string> &variables,
                                                                            const Parameters &at) const
{
  std::vector<arma::mat> result(variables.size());
  for (size_t k = 0; k < variables.size(); ++k)
  {
    const arma::mat &current_parameter = at[variables[k]];
    result[k].zeros(current_parameter.n_rows, current_parameter.n_cols);
  }
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      std::vector<arma::mat> factor_gradients = this->likelihood_estimator_outputs[*i]->get_gradients_of_log(variables, at);
      for (size_t k = 0; k < variables.size(); ++k)
        result[k] = result[k] + factor_gradients[k];
    }
  }
  return result;
}

// (C222) each factor in one pass where it can, else evaluated and asked for
// its gradients separately; the values summed in the index's order, as
// evaluate_likelihoods
double VectorFactorVariables::evaluate_likelihoods_with_gradients(const Index* index,
                                                                  const std::vector<std::string> &variables,
                                                                  std::vector<arma::mat> &gradients)
{
  const Parameters &at = this->particle->parameters;
  gradients.resize(variables.size());
  for (size_t k = 0; k < variables.size(); ++k)
  {
    const arma::mat &current_parameter = at[variables[k]];
    gradients[k].zeros(current_parameter.n_rows, current_parameter.n_cols);
  }
  double result = 0.0;
  std::vector<arma::mat> factor_gradients;
  for (auto i = index->begin(); i != index->end(); ++i)
  {
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    if (output == NULL)
      continue;
    double value;
    if (!output->evaluate_with_gradients(at, variables, value, factor_gradients))
    {
      value = output->evaluate(at);
      factor_gradients = output->get_gradients_of_log(variables, at);
    }
    result = result + value;
    for (size_t k = 0; k < variables.size(); ++k)
      gradients[k] = gradients[k] + factor_gradients[k];
  }
  return result;
}

arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                      const std::string &variable) const
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      result = result + this->likelihood_estimator_outputs[*i]->subsample_get_gradient_of_log(variable,
                                                                                              this->particle->parameters);
    }
  }
  return result;
}

const Factors* VectorFactorVariables::get_factors() const
{
  return this->vector_factors;
}

void VectorFactorVariables::forget_you_were_already_written_to_file()
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    this->likelihood_estimator_outputs[i]->forget_you_were_already_written_to_file();
  }
}

bool VectorFactorVariables::recycle_all_outputs(const Factors* factors,
                                                const Parameters &new_parameters)
{
  // what factors->simulate_factor_variables(new_parameters) returns: every
  // output initialised and simulated at new_parameters. All outputs are
  // reinitialised before any is simulated, so that a failure (an output that
  // cannot be reused) happens before any random numbers are used: the caller
  // then makes new objects, which use exactly the random numbers new objects
  // would.
  if ((factors == NULL) || (factors != this->vector_factors))
    return false;
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[i].get();
    if ((output == NULL) || !output->reinitialise(new_parameters))
      return false;
  }
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    this->likelihood_estimator_outputs[i]->simulate(new_parameters);
    this->likelihood_estimator_outputs[i]->write_to_file_flag = false;
  }
  return true;
}

bool VectorFactorVariables::recycle_with_new_outputs(const FactorVariables* source_in,
                                                     const Parameters &new_parameters,
                                                     const Index* new_outputs_index)
{
  const VectorFactorVariables* source = dynamic_cast<const VectorFactorVariables*>(source_in);
  if ((source == NULL) || (source->vector_factors != this->vector_factors) ||
      (source->likelihood_estimator_outputs.size() != this->likelihood_estimator_outputs.size()))
    return false;
  // as copy_with_new_outputs, output by output
  std::vector<std::string> changed;
  int changed_known = -1; // (found when first needed)
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[i].get();
    bool new_output = std::find(new_outputs_index->begin(),
                                new_outputs_index->end(),
                                i)!=new_outputs_index->end();
    if (new_output)
    {
      if ((output == NULL) || !output->reinitialise(new_parameters))
        return false;
      output->simulate(new_parameters);
      if ((source->likelihood_estimator_outputs[i]!=NULL) && (source->particle!=NULL) && output->carries_over())
      {
        if (changed_known < 0)
          changed_known = changed_variables(source->particle->parameters, new_parameters, changed) ? 1 : 0;
        output->carry_over_from(*source->likelihood_estimator_outputs[i],
                                source->particle->parameters,
                                new_parameters,
                                (changed_known == 1) ? &changed : NULL);
      }
    }
    else if (source->likelihood_estimator_outputs[i]!=NULL)
    {
      if ((output == NULL) || !output->assign_from(*source->likelihood_estimator_outputs[i]))
        return false;
    }
    else if (output != NULL)
    {
      return false;
    }
  }
  return true;
}

int VectorFactorVariables::adaptive_summary_size(const Index* index) const
{
  if (index==NULL)
    return -1;
  int size = 0;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    const LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    if (output==NULL)
      return -1;
    int n = output->adaptive_summary_size();
    if (n < 0)
      return -1;
    size = size + n;
  }
  return size;
}

void VectorFactorVariables::adaptive_summary(const Index* index, double* values) const
{
  size_t offset = 0;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    const LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    output->adaptive_summary(values + offset);
    offset = offset + size_t(output->adaptive_summary_size());
  }
}

double VectorFactorVariables::smcadaptive_from_summary(const Index* index,
                                                       const double* values,
                                                       const Parameters &parameters) const
{
  // as evaluate_smcadaptive_part_given_smcfixed_likelihoods: the same sum in
  // the same order
  double result = 0.0;
  size_t offset = 0;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    const LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    result = result + output->smcadaptive_from_summary(values + offset, parameters);
    offset = offset + size_t(output->adaptive_summary_size());
  }
  return result;
}

bool VectorFactorVariables::restore_from_summary(const Index* index,
                                                const Index* other_index,
                                                const double* values,
                                                const Parameters &parameters)
{
  size_t offset = 0;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    if (output==NULL)
      return false;
    output->restore_from_summary(values + offset, parameters);
    offset = offset + size_t(output->adaptive_summary_size());
  }
  if (other_index!=NULL)
  {
    for (auto i=other_index->begin(); i!=other_index->end(); ++i)
    {
      LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
      if ((output==NULL) || !output->reinitialise(this->particle->parameters))
        return false;
      output->simulate(this->particle->parameters);
      output->write_to_file_flag = false;
    }
  }
  return true;
}

bool VectorFactorVariables::repoint_to(const Factors* factors)
{
  const VectorFactors* vector_factors_in = dynamic_cast<const VectorFactors*>(factors);
  if ((vector_factors_in == NULL) ||
      (vector_factors_in->number_of_likelihood_estimators() != this->likelihood_estimator_outputs.size()))
    return false;
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    if ((this->likelihood_estimator_outputs[i] != NULL) &&
        !this->likelihood_estimator_outputs[i]->repoint_estimator(vector_factors_in->get_likelihood_estimator(i)))
      return false;
  }
  this->vector_factors = vector_factors_in;
  return true;
}

void VectorFactorVariables::smcadaptive_from_summaries(const Index* index,
                                                         const double* values,
                                                         size_t stride,
                                                         size_t n,
                                                         const Parameters &parameters,
                                                         double* results) const
{
  // as smcadaptive_from_summary for each point: the same sum, in the same
  // order (results start at 0.0 and each output's value is added in turn)
  std::fill(results, results + n, 0.0);
  std::vector<double> values_of_output(n);
  size_t offset = 0;
  for (auto i=index->begin(); i!=index->end(); ++i)
  {
    const LikelihoodEstimatorOutput* output = this->likelihood_estimator_outputs[*i].get();
    output->smcadaptive_from_summaries(values + offset, stride, n, parameters, values_of_output.data());
    for (size_t j = 0; j < n; ++j)
      results[j] = results[j] + values_of_output[j];
    offset = offset + size_t(output->adaptive_summary_size());
  }
}

bool VectorFactorVariables::has_compact_values() const
{
  for (auto i = this->likelihood_estimator_outputs.begin(); i != this->likelihood_estimator_outputs.end(); ++i)
  {
    if ((*i)->compact_size() == 0)
      return false;
  }
  return true;
}

size_t VectorFactorVariables::number_of_compact_values() const
{
  size_t n = 0;
  for (auto i = this->likelihood_estimator_outputs.begin(); i != this->likelihood_estimator_outputs.end(); ++i)
    n += (*i)->compact_size();
  return n;
}

void VectorFactorVariables::compact_values(double *values) const
{
  size_t offset = 0;
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    this->likelihood_estimator_outputs[i]->compact_values(values + offset);
    offset += this->likelihood_estimator_outputs[i]->compact_size();
  }
}

void VectorFactorVariables::write_compact_values(const std::string &directory_name,
                                                 const std::string &index,
                                                 const double *values) const
{
  // the same directories as write_to_file, built once per directory
  if (!this->factor_directories)
    this->factor_directories.reset(new FactorDirectories());
  FactorDirectories &cache = *this->factor_directories;
  if ( (cache.names.size() != this->likelihood_estimator_outputs.size()) ||
       (cache.directory != directory_name) )
  {
    cache.directory = directory_name;
    cache.names.clear();
    for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
      cache.names.push_back(directory_name + "/factor" + std::to_string(i+1));
  }
  size_t offset = 0;
  for (size_t i=0; i<this->likelihood_estimator_outputs.size(); ++i)
  {
    this->likelihood_estimator_outputs[i]->write_compact_values(cache.names[i], index, values + offset);
    offset += this->likelihood_estimator_outputs[i]->compact_size();
  }
}

void VectorFactorVariables::write_to_file(const std::string &directory_name,
                                          const std::string &index) const
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    std::string factor_directory_name = directory_name + "/factor" + std::to_string(i+1);
    this->likelihood_estimator_outputs[i]->write_to_file(factor_directory_name,
                                                         index);
  }
}

void VectorFactorVariables::close_ofstreams()
{
  for (size_t i=0;
       i<this->likelihood_estimator_outputs.size();
       ++i)
  {
    this->likelihood_estimator_outputs[i]->close_ofstreams();
  }
}
LikelihoodEstimatorOutput* VectorFactorVariables::get_output(size_t i) const
{
  return this->likelihood_estimator_outputs[i].get();
}

std::unique_ptr<LikelihoodEstimatorOutput> VectorFactorVariables::replace_with_new_output(size_t i,
                                                                                          const Parameters &parameters)
{
  std::unique_ptr<LikelihoodEstimatorOutput> old = std::move(this->likelihood_estimator_outputs[i]);
  LikelihoodEstimatorOutput* output = this->vector_factors->get_likelihood_estimator(i)->initialise(parameters);
  output->simulate(parameters);
  output->write_to_file_flag = false;
  this->likelihood_estimator_outputs[i].reset(output);
  return old;
}

void VectorFactorVariables::put_back_output(size_t i, std::unique_ptr<LikelihoodEstimatorOutput> output)
{
  this->likelihood_estimator_outputs[i] = std::move(output);
}

}
