#include "vector_factor_variables.h"
#include "likelihood_estimator_output.h"
#include "index.h"
#include "particle.h"
#include "vector_factors.h"

VectorFactorVariables::VectorFactorVariables()
  :FactorVariables()
{
  this->likelihood_estimator_outputs.resize(0);
}

VectorFactorVariables::~VectorFactorVariables()
{
  for (std::vector<LikelihoodEstimatorOutput*>::iterator i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

VectorFactorVariables::VectorFactorVariables(VectorFactors* vector_factors_in,
                                             const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in)
:FactorVariables()
{
  this->vector_factors = vector_factors_in;
  this->likelihood_estimator_outputs = likelihood_estimator_outputs_in;
}

VectorFactorVariables::VectorFactorVariables(VectorFactors* vector_factors_in,
                                             const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                                             Particle* particle_in)
:FactorVariables(particle_in)
{
  this->vector_factors = vector_factors_in;
  this->likelihood_estimator_outputs = likelihood_estimator_outputs_in;
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
  
  for (std::vector<LikelihoodEstimatorOutput*>::iterator i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimator_outputs.clear();
  
  FactorVariables::operator=(another);
  this->make_copy(another);
}

FactorVariables* VectorFactorVariables::duplicate() const
{
  return( new VectorFactorVariables(*this));
}

void VectorFactorVariables::make_copy(const VectorFactorVariables &another)
{
  this->likelihood_estimator_outputs.resize(0);
  this->likelihood_estimator_outputs.reserve(another.likelihood_estimator_outputs.size());
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=another.likelihood_estimator_outputs.begin();
       i!=another.likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimator_outputs.push_back((*i)->duplicate());
    else
      this->likelihood_estimator_outputs.push_back(NULL);
  }
  this->vector_factors = another.vector_factors;
}

void VectorFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index)
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

void VectorFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                                  const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(all_parameters);
    }
  }
}

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

double VectorFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    if (this->likelihood_estimator_outputs[*i]!=NULL)
    {
      this->likelihood_estimator_outputs[*i]->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
      result = result + this->likelihood_estimator_outputs[*i]->log_likelihood;
    }
  }
  //this->target_evaluated = result;
  return result;
}

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

double VectorFactorVariables::evaluate_likelihoods(const Index* index)
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

double VectorFactorVariables::subsample_evaluate_likelihoods(const Index* index)
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

arma::mat VectorFactorVariables::direct_get_gradient_of_log(const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=this->likelihood_estimator_outputs.begin();
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

arma::mat VectorFactorVariables::direct_get_gradient_of_log(const std::string &variable,
                                               const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=this->likelihood_estimator_outputs.begin();
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

arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable)
{
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=this->likelihood_estimator_outputs.begin();
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

arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                                      const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->particle->parameters[variable];
  arma::mat result(current_parameter.n_rows,current_parameter.n_cols);
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=this->likelihood_estimator_outputs.begin();
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

arma::mat VectorFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                            const std::string &variable)
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

arma::mat VectorFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                      const std::string &variable)
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

Factors* VectorFactorVariables::get_factors()
{
  return this->vector_factors;
}
