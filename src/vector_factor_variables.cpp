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

FactorVariables* VectorFactorVariables::copy_with_new_outputs(const Parameters &new_parameters,
                                                             const Index* new_outputs_index) const
{
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimator_outputs.size());
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
    //std::cout << "Index is NULL." << std::endl;
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

/*
 void VectorFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters;
 if (conditioned_on_parameters.is_empty())
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
 else if (this->particle->parameters.is_empty())
 {
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(conditioned_on_parameters);
 }
 }
 }
 else
 {
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(this->particle->parameters.merge(conditioned_on_parameters));
 }
 }
 }
 }
 */

void VectorFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  //Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
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

/*
 void VectorFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcfixed_part(all_parameters);
 }
 }
 }
 */

double VectorFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    //if (this->likelihood_estimator_outputs[*i]!=NULL)
    //{
    
    this->likelihood_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
    result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
    
    //}
  }
  //this->target_evaluated = result;
  return result;
}

/*
 double VectorFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 double result = 0.0;
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 this->likelihood_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(all_parameters);
 result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
 }
 }
 //this->target_evaluated = result;
 return result;
 }
 */

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
  //this->subsample_target_evaluated = result;
  return result;
}

/*
 double VectorFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 double result = 0.0;
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 this->likelihood_estimator_outputs[*i]->subsample_evaluate_smcadaptive_part_given_smcfixed(all_parameters);
 result = result + this->likelihood_estimator_outputs[*i]->subsample_log_likelihood;
 }
 }
 //this->subsample_target_evaluated = result;
 return result;
 }
 */

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
  //this->target_evaluated = result;
  return result;
}

/*
 double VectorFactorVariables::evaluate_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 double result = 0.0;
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 result = result + this->likelihood_estimator_outputs[*i]->evaluate(all_parameters);
 }
 }
 //this->target_evaluated = result;
 return result;
 }
 */

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
  //this->subsample_target_evaluated = result;
  return result;
}

/*
 double VectorFactorVariables::subsample_evaluate_likelihoods(const Index* index,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 double result = 0.0;
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 result = result + this->likelihood_estimator_outputs[*i]->subsample_evaluate(all_parameters);
 }
 }
 //this->subsample_target_evaluated = result;
 return result;
 }
 */

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
  //this->target_gradients_of_log[variable] = result;
  return result;
}

/*
 arma::mat VectorFactorVariables::direct_get_gradient_of_log(const std::string &variable,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 arma::mat current_parameter = this->particle->parameters[variable];
 arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
 for (auto i=this->likelihood_estimator_outputs.begin();
 i!=this->likelihood_estimator_outputs.end();
 ++i)
 {
 if (*i!=NULL)
 {
 result = result + (*i)->get_gradient_of_log(variable,
 all_parameters);
 }
 }
 //this->target_gradients_of_log[variable] = result;
 return result;
 }
 */

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
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

/*
 arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 arma::mat current_parameter = this->particle->parameters[variable];
 arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
 for (auto i=this->likelihood_estimator_outputs.begin();
 i!=this->likelihood_estimator_outputs.end();
 ++i)
 {
 if (*i!=NULL)
 {
 result = result + (*i)->subsample_get_gradient_of_log(variable,
 all_parameters);
 }
 }
 //this->subsample_target_gradients_of_log[variable] = result;
 return result;
 }
 */

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
  //this->target_gradients_of_log[variable] = result;
  return result;
}

/*
 arma::mat VectorFactorVariables::direct_get_gradient_of_log(const Index* index,
 const std::string &variable,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 arma::mat current_parameter = this->particle->parameters[variable];
 arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 result = result + this->likelihood_estimator_outputs[*i]->get_gradient_of_log(variable,
 all_parameters);
 }
 }
 //this->target_gradients_of_log[variable] = result;
 return result;
 }
 */

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
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

/*
 arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
 const std::string &variable,
 const Parameters &conditioned_on_parameters)
 {
 Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
 arma::mat current_parameter = this->particle->parameters[variable];
 arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
 for (auto i=index->begin();
 i!=index->end();
 ++i)
 {
 if (this->likelihood_estimator_outputs[*i]!=NULL)
 {
 result = result + this->likelihood_estimator_outputs[*i]->subsample_get_gradient_of_log(variable,
 all_parameters);
 }
 }
 //this->subsample_target_gradients_of_log[variable] = result;
 return result;
 }
 */

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
}
