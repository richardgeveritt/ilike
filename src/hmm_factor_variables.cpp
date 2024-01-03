#include "hmm_factor_variables.h"
#include "hmm_factors.h"
#include "index.h"
#include "likelihood_estimator_output.h"
#include "particle.h"
#include "proposal_kernel.h"

HMMFactorVariables::HMMFactorVariables()
  :FactorVariables()
{
  this->likelihood_estimator_outputs.resize(0);
  //this->initial_prior = NULL;
  this->dynamic_smcfixed_part = 0.0;
  this->hmm_factors = NULL;
}

HMMFactorVariables::HMMFactorVariables(const HMMFactors* hmm_factors_in,
                                       const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in)
:FactorVariables()
{
  this->hmm_factors = hmm_factors_in;
  this->likelihood_estimator_outputs = likelihood_estimator_outputs_in;
}

HMMFactorVariables::HMMFactorVariables(const HMMFactors* hmm_factors_in,
                                       const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                                       Particle* particle_in)
:FactorVariables(particle_in)
{
  this->hmm_factors = hmm_factors_in;
  this->likelihood_estimator_outputs = likelihood_estimator_outputs_in;
}

HMMFactorVariables::~HMMFactorVariables()
{
  for (std::vector<LikelihoodEstimatorOutput*>::iterator i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  //if (this->initial_prior!=NULL)
  //  delete this->initial_prior;
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
  
  for (std::vector<LikelihoodEstimatorOutput*>::iterator i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->likelihood_estimator_outputs.clear();
  
  //if (this->initial_prior!=NULL)
  //  delete this->initial_prior;
  
  FactorVariables::operator=(another);
  this->make_copy(another);
}

FactorVariables* HMMFactorVariables::duplicate() const
{
  return( new HMMFactorVariables(*this));
}

void HMMFactorVariables::make_copy(const HMMFactorVariables &another)
{
  this->hmm_factors = another.hmm_factors;
  
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
  
  /*
  if (another.initial_prior!=NULL)
    this->initial_prior = another.initial_prior->duplicate();
  else
    this->initial_prior = NULL;
  */
}

void HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  //if (index->size()>1)
  //  Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  /*
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
      this->initial_prior->evaluate_smcfixed_part(this->particle->parameters);
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL))
      this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                                          *this->particle->previous_self);
  }
  */
  
  if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL) && (index->get_transition_model()))
    this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                                        *this->particle->previous_self);
  
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    this->likelihood_estimator_outputs[*i]->evaluate_smcfixed_part(this->particle->parameters);
  }
  
  //if (this->likelihood_estimator_outputs[*index.begin()]!=NULL)
  //{
  //  this->likelihood_estimator_outputs->
  //  this->likelihood_estimator_outputs[*index.begin()]->evaluate_smcfixed_part(this->particle->parameters);
  //}
}

/*
void HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
      this->initial_prior->evaluate_smcfixed_part(all_parameters);
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL))
      this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                                          *this->particle->previous_self,
                                                                                          conditioned_on_parameters);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->evaluate_smcfixed_part(all_parameters);
  }
}
*/

void HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  
  /*
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
      this->initial_prior->evaluate_smcfixed_part(this->particle->parameters);
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL))
      this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                                    *this->particle->previous_self);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->subsample_evaluate_smcfixed_part(this->particle->parameters);
  }
  */
  
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

/*
void HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                                         const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
      this->initial_prior->evaluate_smcfixed_part(all_parameters);
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag && (this->hmm_factors->transition_kernel!=NULL))
      this->dynamic_smcfixed_part = this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                                    *this->particle->previous_self,
                                                                                                    conditioned_on_parameters);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->subsample_evaluate_smcfixed_part(all_parameters);
  }
}
*/

double HMMFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  //if (index->size()>1)
  //  Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  double result = 0.0;
  
  /*
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      this->initial_prior->evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
      result = result + this->initial_prior->log_likelihood;
    }
  }
  else
  {
    
  }
  */
  
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

/*
double HMMFactorVariables::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  double result = 0.0;
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      this->initial_prior->evaluate_smcadaptive_part_given_smcfixed(all_parameters);
      result = result + this->initial_prior->log_likelihood;
    }
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag)
      result = result + this->dynamic_smcfixed_part;
    else
    {
      if (this->hmm_factors->transition_kernel!=NULL)
        result = result + this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                                *this->particle->previous_self,
                                                                                conditioned_on_parameters);
    }
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->evaluate_smcadaptive_part_given_smcfixed(all_parameters);
    result = result + (*i)->log_likelihood;
  }
  
  return result;
}
*/

/*
double HMMFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                          const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  double result = 0.0;
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      this->initial_prior->subsample_evaluate_smcadaptive_part_given_smcfixed(all_parameters);
      result = result + this->initial_prior->subsample_log_likelihood;
    }
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag)
      result = result + this->dynamic_smcfixed_part;
    else
    {
      if (this->hmm_factors->transition_kernel!=NULL)
        result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                          *this->particle->previous_self,
                                                                                          conditioned_on_parameters);
    }
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->subsample_evaluate_smcadaptive_part_given_smcfixed(all_parameters);
    result = result + (*i)->subsample_log_likelihood;
  }
  
  return result;
}
*/

double HMMFactorVariables::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  
  /*
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_smcfixed_part_of_likelihoods - cannot break down likelihood estimation into two stages when using more than one measurement.");
  
  double result = 0.0;
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      this->initial_prior->subsample_evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
      result = result + this->initial_prior->subsample_log_likelihood;
    }
  }
  else
  {
    if (this->hmm_factors->smcfixed_flag)
      result = result + this->dynamic_smcfixed_part;
    else
    {
      if (this->hmm_factors->transition_kernel!=NULL)
        result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                          *this->particle->previous_self);
    }
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    (*i)->subsample_evaluate_smcadaptive_part_given_smcfixed(this->particle->parameters);
    result = result + (*i)->subsample_log_likelihood;
  }
  */
  
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

double HMMFactorVariables::evaluate_likelihoods(const Index* index) const
{
  /*
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::evaluate_likelihoods - only one index allowed.");
  
  double result = 0.0;
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->evaluate(this->particle->parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                              *this->particle->previous_self);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    result = result + (*i)->evaluate(this->particle->parameters);
  }

  return result;
  */
  
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

/*
double HMMFactorVariables::evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::evaluate_likelihoods - only one index allowed.");
  
  double result = 0.0;
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->evaluate(all_parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->evaluate_kernel(*this->particle,
                                                                              *this->particle->previous_self,
                                                                              conditioned_on_parameters);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    result = result + (*i)->evaluate(all_parameters);
  }
  
  return result;
}
*/

double HMMFactorVariables::subsample_evaluate_likelihoods(const Index* index) const
{
  /*
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - only one index allowed.");
  
  double result = 0.0;
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->subsample_evaluate(this->particle->parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                        *this->particle->previous_self);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    result = result + (*i)->subsample_evaluate(this->particle->parameters);
  }
  
  return result;
  */
  
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

/*
double HMMFactorVariables::subsample_evaluate_likelihoods(const Index* index,
                                                          const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->particle->parameters.merge(conditioned_on_parameters);
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::evaluate_likelihoods - only one index allowed.");
  
  double result = 0.0;
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->subsample_evaluate(all_parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->subsample_evaluate_kernel(*this->particle,
                                                                                        *this->particle->previous_self,
                                                                                        conditioned_on_parameters);
  }
  
  for (auto i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    result = result + (*i)->subsample_evaluate(all_parameters);
  }
  
  return result;
}
*/

/*
arma::mat HMMFactorVariables::direct_get_gradient_of_log(const std::string &variable)
{
  if (index->size()>1)
    Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - not written yet since unrolled HMM not yet implemented.");
  
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
  
  if (index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->subsample_evaluate(all_parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->gradient_of_log(variable,
                                                                              *this->particle,
                                                                              *this->particle->previous_self);
  }
  
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMFactorVariables::direct_get_gradient_of_log(const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - not written yet since unrolled HMM not yet implemented.");
  
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
  
  
  if (this->hmm_factors->transition_kernel!=NULL)
    result = result + this->hmm_factors->transition_kernel->gradient_of_log(variable,
                                                                            *this->particle,
                                                                            *this->particle->previous_self,
                                                                            conditioned_on_parameters);
  
  //this->target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable)
{
  Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - not written yet since unrolled HMM not yet implemented.");
  
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
  
  if (this->hmm_factors->transition_kernel!=NULL)
    result = result + this->hmm_factors->transition_kernel->subsample_gradient_of_log(variable,
                                                                                      *this->particle,
                                                                                      *this->particle->previous_self);
  
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}

arma::mat HMMFactorVariables::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                                   const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - not written yet since unrolled HMM not yet implemented.");
  
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
  
  if (this->hmm_factors->transition_kernel!=NULL)
    result = result + this->hmm_factors->transition_kernel->subsample_gradient_of_log(variable,
                                                                                      *this->particle,
                                                                                      *this->particle->previous_self,
                                                                                      conditioned_on_parameters);
  
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
}
*/

arma::mat HMMFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                         const std::string &variable) const
{
  Rcpp::stop("HMMFactorVariables::direct_get_gradient_of_log - not written yet..");
  
  /*
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
  
  if (*index->begin()==0)
  {
    if (this->initial_prior!=NULL)
    {
      result = result + this->initial_prior->subsample_evaluate(this->particle->parameters);
    }
  }
  else
  {
    if (this->hmm_factors->transition_kernel!=NULL)
      result = result + this->hmm_factors->transition_kernel->gradient_of_log(variable,
                                                                              *this->particle,
                                                                              *this->particle->previous_self);
  }
  
  //this->target_gradients_of_log[variable] = result;
  return result;
  */
}

/*
arma::mat HMMFactorVariables::direct_get_gradient_of_log(const Index* index,
                                                         const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("HMMFactorVariables::direct_get_gradient_of_log - not written yet.");
  
}
*/

arma::mat HMMFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                   const std::string &variable) const
{
  Rcpp::stop("HMMFactorVariables::direct_get_gradient_of_log - not written yet.");
  
  /*
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
  
  if (this->hmm_factors->transition_kernel!=NULL)
    result = result + this->hmm_factors->transition_kernel->subsample_gradient_of_log(variable,
                                                                                      *this->particle,
                                                                                      *this->particle->previous_self);
  
  //this->subsample_target_gradients_of_log[variable] = result;
  return result;
  */
}

/*
arma::mat HMMFactorVariables::direct_subsample_get_gradient_of_log(const Index* index,
                                                                   const std::string &variable,
                                                                   const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("HMMFactorVariables::subsample_evaluate_likelihoods - not written yet.");
}
*/

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
