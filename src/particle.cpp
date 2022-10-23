#include "particle.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"
#include "transform.h"
#include "proposal_kernel.h"
#include "gradient_estimator_output.h"
#include "gradient_estimator.h"
#include "factor_variables.h"
#include "vector_factor_variables.h"
#include "ensemble_factor_variables.h"
#include "factors.h"
#include "ensemble_factors.h"

Particle::Particle()
{
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
  
  this->factor_variables = NULL;
  this->ensemble_factor_variables = NULL;
  
  this->move_transformed_parameters = Parameters();
  this->move_transform = NULL;
  this->move_parameters = NULL;
}

Particle::Particle(const Parameters &parameters_in)
{
  this->parameters = parameters_in;
  this->move_transformed_parameters = Parameters();
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
  
  this->factor_variables = NULL;
  this->ensemble_factor_variables = NULL;
}

Particle::Particle(const Parameters &parameters_in,
                   FactorVariables* factor_variables_in)
{
  this->parameters = parameters_in;
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  this->factor_variables = factor_variables_in;
  this->ensemble_factor_variables = NULL;
  //std::cout<<this->parameters<<std::endl;
  this->factor_variables->set_particle(this);
  //std::cout<<this->factor_variables->get_particle()->parameters<<std::endl;
}


Particle::Particle(const Parameters &parameters_in,
                   EnsembleFactorVariables* ensemble_factor_variables_in)
{
  this->parameters = parameters_in;
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  this->ensemble_factor_variables = ensemble_factor_variables_in;
  this->factor_variables = NULL;
}

Particle::Particle(const Parameters &parameters_in,
                   FactorVariables* factor_variables_in,
                   double previous_target_evaluated_in)
{
  this->parameters = parameters_in;
  this->previous_target_evaluated = previous_target_evaluated_in;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
  
  this->factor_variables = factor_variables_in;
  this->factor_variables->set_particle(this);
  this->ensemble_factor_variables = NULL;
}

Particle::~Particle()
{
  if (this->factor_variables!=NULL)
    delete this->factor_variables;
  
  if (this->ensemble_factor_variables!=NULL)
    delete this->ensemble_factor_variables;
  
  for (auto i=this->gradient_estimator_outputs.begin();
       i!=this->gradient_estimator_outputs.end();
       ++i)
  {
    if (i->second!=NULL)
      delete i->second;
  }
}

//Copy constructor for the Particle class.
Particle::Particle(const Particle &another)
{
  this->make_copy(another);
}

void Particle::operator=(const Particle &another)
{
  if(this == &another)
    return;
  
  if (this->factor_variables!=NULL)
    delete this->factor_variables;
  
  if (this->ensemble_factor_variables!=NULL)
    delete this->ensemble_factor_variables;
  
  for (auto i=this->gradient_estimator_outputs.begin();
       i!=this->gradient_estimator_outputs.end();
       ++i)
  {
    if (i->second!=NULL)
      delete i->second;
  }
  this->gradient_estimator_outputs.clear();

  this->make_copy(another);
}

void Particle::make_copy(const Particle &another)
{
  if (another.factor_variables!=NULL)
  {
    this->factor_variables = another.factor_variables->duplicate();
    this->factor_variables->set_particle(this);
  }
  else
    this->factor_variables = NULL;
  
  if (another.ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables = another.ensemble_factor_variables->duplicate();
    this->ensemble_factor_variables->set_particle(this);
  }
  else
    this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = another.previous_target_evaluated;
  this->target_evaluated = another.target_evaluated;
  this->subsample_previous_target_evaluated = another.subsample_previous_target_evaluated;
  this->subsample_target_evaluated = another.subsample_target_evaluated;
  
  this->previous_ensemble_target_evaluated = another.previous_ensemble_target_evaluated;
  this->subsample_previous_ensemble_target_evaluated = another.subsample_previous_ensemble_target_evaluated;
  this->ensemble_target_evaluated = another.ensemble_target_evaluated;
  this->subsample_ensemble_target_evaluated = another.subsample_ensemble_target_evaluated;
  
  this->parameters = another.parameters;
  
  this->gradient_estimator_outputs = boost::unordered_map<const ProposalKernel*, GradientEstimatorOutput*>();
  this->gradient_estimator_outputs.reserve(another.gradient_estimator_outputs.size());
  for (auto i=another.gradient_estimator_outputs.begin();
       i!=another.gradient_estimator_outputs.end();
       ++i)
  {
    if (i->second!=NULL)
      this->gradient_estimator_outputs[i->first] = i->second->duplicate();
    else
      this->gradient_estimator_outputs[i->first] = NULL;
  }
  
  this->move_transformed_parameters = another.move_transformed_parameters;
  this->move_transform = another.move_transform;
  this->move_parameters = another.move_parameters;
  this->accepted_outputs = another.accepted_outputs;
  
  this->previous_self = another.previous_self;
}

/*
void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index);
}

void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                     const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index,
                                                                  conditioned_on_parameters);
}

void Particle::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->subsample_evaluate_smcfixed_part_of_likelihoods(index,
                                                                            conditioned_on_parameters);
}

double Particle::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                      const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                                          conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                                                    conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_likelihoods(index,
                                                                          conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->factor_variables->subsample_evaluate_likelihoods(index);
    return this->subsample_target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->factor_variables->subsample_evaluate_likelihoods(index,
                                                                                              conditioned_on_parameters);
    return this->subsample_target_evaluated;
  }
  else
    return 0.0;
}
*/

void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index);
}

void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                     const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index,
                                                                  conditioned_on_parameters);
}

void Particle::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->subsample_evaluate_smcfixed_part_of_likelihoods(index,
                                                                            conditioned_on_parameters);
}

double Particle::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                      const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                        conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                        conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    //std::cout<<this->parameters<<std::endl;
    //std::cout<<this->factor_variables->get_particle()->parameters<<std::endl;
    this->target_evaluated = this->factor_variables->evaluate_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_likelihoods(index,
                                                        conditioned_on_parameters);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->factor_variables->subsample_evaluate_likelihoods(index);
    return this->subsample_target_evaluated;
  }
  else
    return 0.0;
}

double Particle::subsample_evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->factor_variables->subsample_evaluate_likelihoods(index,
                                                                  conditioned_on_parameters);
    return this->subsample_target_evaluated;
  }
  else
    return 0.0;
}

/*
std::vector<LikelihoodEstimatorOutput*> Particle::get_likelihood_estimator_outputs() const
{
  return this->likelihood_estimator_outputs;
}
*/

arma::mat Particle::direct_get_gradient_of_log(const std::string &variable,
                                               const Index* index)
{
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_get_gradient_of_log(index,variable);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

arma::mat Particle::direct_get_gradient_of_log(const std::string &variable,
                                               const Index* index,
                                               const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_get_gradient_of_log(index,
                                                              variable,
                                                              all_parameters);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

arma::mat Particle::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                         const Index* index)
{
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_subsample_get_gradient_of_log(index,
                                                                        variable);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

arma::mat Particle::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                         const Index* index,
                                                         const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_subsample_get_gradient_of_log(index,
                                                                        variable,
                                                                        all_parameters);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

double Particle::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double incremental_temperature)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                incremental_temperature);
  else
    return 0.0;
}

double Particle::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double incremental_temperature,
                                                     const Parameters &conditioned_on_parameters)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                incremental_temperature,
                                                                                conditioned_on_parameters);
  else
    return 0.0;
}

double Particle::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double incremental_temperature,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                incremental_temperature,
                                                                                conditioned_on_parameters);
  else
    return 0.0;
  
}

double Particle::evaluate_all_likelihoods(const Index* index)
{
  this->target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->target_evaluated + this->factor_variables->evaluate_likelihoods(index);
  }
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_target_evaluated = this->ensemble_target_evaluated + this->ensemble_factor_variables->evaluate_likelihoods(index);
  }
  return this->target_evaluated + this->ensemble_factor_variables->ensemble_factors->get_temperature()*this->ensemble_target_evaluated;
}

double Particle::evaluate_all_likelihoods(const Index* index,
                                          const Parameters &conditioned_on_parameters)
{
  this->target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->target_evaluated + this->factor_variables->evaluate_likelihoods(index,
                                                                                                   conditioned_on_parameters);
  }
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_target_evaluated = this->ensemble_target_evaluated + this->ensemble_factor_variables->evaluate_likelihoods(index,
                                                                                                                              conditioned_on_parameters);
  }
  return this->target_evaluated + this->ensemble_factor_variables->ensemble_factors->get_temperature()*this->ensemble_target_evaluated;
}

double Particle::subsample_evaluate_all_likelihoods(const Index* index,
                                                    const Parameters &conditioned_on_parameters)
{
  this->subsample_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->subsample_target_evaluated + this->factor_variables->subsample_evaluate_likelihoods(index,
                                                                                                                                 conditioned_on_parameters);
  }
  if (this->ensemble_factor_variables!=NULL)
  {
    this->subsample_ensemble_target_evaluated = this->subsample_ensemble_target_evaluated + this->ensemble_factor_variables->subsample_evaluate_likelihoods(index,
                                                                                                                                                            conditioned_on_parameters);
  }
  return this->subsample_target_evaluated + this->ensemble_factor_variables->ensemble_factors->get_temperature()*this->subsample_ensemble_target_evaluated;
}

arma::colvec Particle::get_vector() const
{
  return this->parameters.get_vector();
}

arma::colvec Particle::get_vector(const std::vector<std::string> &state_names) const
{
  return this->parameters.get_vector();
}

GradientEstimatorOutput* Particle::initialise_gradient_estimator_output(const ProposalKernel* proposal,
                                                                         GradientEstimator* gradient_estimator)
{
  auto found = this->gradient_estimator_outputs.find(proposal);
  
  if (found != this->gradient_estimator_outputs.end())
  {
    return found->second;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    GradientEstimatorOutput* current_output = gradient_estimator->initialise();
    this->gradient_estimator_outputs[proposal] = current_output;
    return current_output;
  }
}

void Particle::set_move_transformed_parameters()
{
  this->move_parameters = &this->parameters;
}

void Particle::set_move_transformed_parameters(Transform* transform_in)
{
  // if the transform is not set, or is a different one to used previously
  if ( (this->move_transform==NULL) || (this->move_transform!=transform_in) )
  {
    this->move_transform = transform_in;
    this->move_transformed_parameters = this->move_transform->transform(this->parameters);
  }
  
  this->move_parameters = &this->move_transformed_parameters;
}

void Particle::set_acceptance(const ProposalKernel* proposal_in,
                               bool accepted_in)
{
  this->accepted_outputs[proposal_in] = accepted_in;
}

void Particle::erase_mcmc_adaptation_info()
{
  this->accepted_outputs.clear();
}

void Particle::simulate_factor_variables(Particle* previous_particle)
{
  if (previous_particle->factor_variables!=NULL)
  {
    
    FactorVariables* new_factor_variables = previous_particle->factor_variables->get_factors()->simulate_factor_variables(this->parameters);
    if (this->factor_variables!=NULL)
    {
      delete this->factor_variables;
    }
    
    this->factor_variables = new_factor_variables;
  }
  else
    this->factor_variables = NULL;
}

void Particle::simulate_ensemble_factor_variables(Particle* previous_particle)
{
  if (previous_particle->ensemble_factor_variables!=NULL)
  {
    EnsembleFactorVariables* new_ensemble_factor_variables = previous_particle->ensemble_factor_variables->get_ensemble_factors()->simulate_ensemble_factor_variables(this->parameters);
    
    if (this->ensemble_factor_variables!=NULL)
    {
      delete this->ensemble_factor_variables;
    }
    
    this->ensemble_factor_variables = new_ensemble_factor_variables;
  }
  else
    this->ensemble_factor_variables = NULL;
}

void Particle::simulate_factor_variables(Particle* previous_particle,
                                         const Parameters &conditioned_on_parameters)
{
  if (previous_particle->factor_variables!=NULL)
  {
    Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
    
    FactorVariables* new_factor_variables = previous_particle->factor_variables->get_factors()->simulate_factor_variables(all_parameters);
    if (this->factor_variables!=NULL)
    {
      delete this->factor_variables;
    }
    
    this->factor_variables = new_factor_variables;
  }
  else
    this->factor_variables = NULL;
}

void Particle::simulate_ensemble_factor_variables(Particle* previous_particle,
                                                  const Parameters &conditioned_on_parameters)
{
  if (previous_particle->ensemble_factor_variables!=NULL)
  {
    Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
    
    EnsembleFactorVariables* new_ensemble_factor_variables = previous_particle->ensemble_factor_variables->get_ensemble_factors()->simulate_ensemble_factor_variables(all_parameters);
    
    if (this->ensemble_factor_variables!=NULL)
    {
      delete this->ensemble_factor_variables;
    }
    
    this->ensemble_factor_variables = new_ensemble_factor_variables;
  }
  else
    this->ensemble_factor_variables = NULL;
}

std::ostream& operator<<(std::ostream& os, const Particle &p)
{
  os << p.parameters << std::endl;

  os << "{" << std::endl;

  /*
  std::vector<LikelihoodEstimatorOutput*>::const_iterator it;

  for (it=p.likelihood_estimator_outputs.begin();it!=p.likelihood_estimator_outputs.end();++it)
  {
    if (it==p.likelihood_estimator_outputs.begin())
      os << **it << std::endl;
    else
      os << "," << *  *it << std::endl;
  }
   */

  os << "}";

  return os;
}
