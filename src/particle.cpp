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
: parameters(),
move_transformed_parameters()
{
  this->move_transform = NULL;
  this->move_parameters = NULL;
  
  this->factor_variables = NULL;
  this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(const Parameters &parameters_in,
                   EnsembleFactors* ensemble_factors_in)
: parameters(parameters_in),
move_transformed_parameters()
{
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters);
    this->ensemble_factor_variables->set_particle(this);
  }
  
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   Factors* factors_in)
: parameters(std::move(parameters_in)),
move_transformed_parameters()
{
  if (factors_in!=NULL)
  {
    this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
    this->factor_variables->set_particle(this);
  }
  
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;

  this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   Factors* factors_in,
                   const Parameters &conditioned_on_parameters)
: parameters(std::move(parameters_in)),
move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
    this->factor_variables->set_particle(this);
  }
  
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   Factors* factors_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
: parameters(std::move(parameters_in)),
move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  this->parameters.merge_with_fixed(sequencer_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
    this->factor_variables->set_particle(this);
  }
    
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;

  this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   FactorVariables* factor_variables_in)
: parameters(std::move(parameters_in)),
move_transformed_parameters(),
factor_variables(factor_variables_in)
{
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->ensemble_factor_variables = NULL;

  this->factor_variables->set_particle(this);
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   FactorVariables* factor_variables_in,
                   double previous_target_evaluated_in)
: parameters(std::move(parameters_in)),
move_transformed_parameters(),
factor_variables(factor_variables_in)
{
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->ensemble_factor_variables = NULL;
  
  this->factor_variables->set_particle(this);
  
  this->previous_target_evaluated = previous_target_evaluated_in;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   EnsembleFactors* ensemble_factors_in)
: parameters(std::move(parameters_in)),
move_transformed_parameters()
{
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters);
    this->ensemble_factor_variables->set_particle(this);
  }
  
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;

  this->factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters)
{
  this->parameters = std::move(parameters_in);
  
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->move_transformed_parameters = Parameters();
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters);
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables = NULL;
  }
  this->factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
{
  this->parameters = std::move(parameters_in);
  
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  this->parameters.merge_with_fixed(sequencer_parameters);
  
  this->move_transformed_parameters = Parameters();
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters);
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables = NULL;
  }
  this->factor_variables = NULL;
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(Parameters &&parameters_in,
                   EnsembleFactorVariables* ensemble_factor_variables_in)
{
  this->parameters = std::move(parameters_in);
  this->move_transformed_parameters = Parameters();
  this->move_transform = NULL;
  this->move_parameters = &this->parameters;
  
  this->ensemble_factor_variables = ensemble_factor_variables_in;
  this->factor_variables = NULL;
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables->set_particle(this);
  }
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

void Particle::setup(Factors* factors_in)
{
  this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(Factors* factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(Factors* factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  this->parameters.merge_with_fixed(sequencer_parameters);
  
  this->factor_variables = factors_in->simulate_factor_variables(this->parameters);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in)
{
  this->parameters = parameters_in;
  this->factor_variables = factors_in->simulate_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->factor_variables = factors_in->simulate_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  this->parameters.merge_with_fixed(sequencer_parameters);
  
  this->factor_variables = factors_in->simulate_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in)
{
  this->parameters = parameters_in;
  this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->ensemble_factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->ensemble_factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  this->parameters.merge_with_fixed(sequencer_parameters);
  
  this->ensemble_factor_variables = ensemble_factors_in->simulate_ensemble_factor_variables(parameters_in);
  this->move_parameters = &this->parameters;
  this->ensemble_factor_variables->set_particle(this);
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

Particle& Particle::operator=(const Particle &another)
{
  if(this == &another)
    return *this;
  
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
  
  return *this;
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
  
  if (this->move_parameters==&another.parameters)
  {
    this->move_parameters = &this->parameters;
  }
  else if (this->move_parameters==&another.move_transformed_parameters)
  {
    this->move_parameters = &this->move_transformed_parameters;
  }
  
  this->accepted_outputs = another.accepted_outputs;
  
  this->previous_self = another.previous_self;
}

//Move constructor for the Particle class.
Particle::Particle(Particle &&another)
{
  this->make_copy(std::move(another));
}

Particle& Particle::operator=(Particle &&another)
{
  if(this == &another)
    return *this;
  
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
  
  this->make_copy(std::move(another));
  
  return *this;
}

void Particle::make_copy(Particle &&another)
{
  if (another.factor_variables!=NULL)
  {
    this->factor_variables = another.factor_variables;
    this->factor_variables->set_particle(this);
  }
  else
    this->factor_variables = NULL;
  
  if (another.ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables = another.ensemble_factor_variables;
    this->ensemble_factor_variables->set_particle(this);
  }
  else
    this->ensemble_factor_variables = NULL;
  
  this->previous_target_evaluated = std::move(another.previous_target_evaluated);
  this->target_evaluated = std::move(another.target_evaluated);
  this->subsample_previous_target_evaluated = std::move(another.subsample_previous_target_evaluated);
  this->subsample_target_evaluated = std::move(another.subsample_target_evaluated);
  
  this->previous_ensemble_target_evaluated = std::move(another.previous_ensemble_target_evaluated);
  this->subsample_previous_ensemble_target_evaluated = std::move(another.subsample_previous_ensemble_target_evaluated);
  this->ensemble_target_evaluated = std::move(another.ensemble_target_evaluated);
  this->subsample_ensemble_target_evaluated = std::move(another.subsample_ensemble_target_evaluated);
  
  this->parameters = std::move(another.parameters);
  
  this->gradient_estimator_outputs = boost::unordered_map<const ProposalKernel*, GradientEstimatorOutput*>();
  this->gradient_estimator_outputs.reserve(another.gradient_estimator_outputs.size());
  for (auto i=another.gradient_estimator_outputs.begin();
       i!=another.gradient_estimator_outputs.end();
       ++i)
  {
    if (i->second!=NULL)
      this->gradient_estimator_outputs[i->first] = i->second;
    else
      this->gradient_estimator_outputs[i->first] = NULL;
  }
  
  this->move_transformed_parameters = std::move(another.move_transformed_parameters);
  this->move_transform = std::move(another.move_transform);
  
  if (this->move_parameters==&another.parameters)
  {
    this->move_parameters = &this->parameters;
  }
  else if (this->move_parameters==&another.move_transformed_parameters)
  {
    this->move_parameters = &this->move_transformed_parameters;
  }
  
  this->accepted_outputs = std::move(another.accepted_outputs);
  
  this->previous_self = another.previous_self;
  
  another.factor_variables = NULL;
  another.ensemble_factor_variables = NULL;
  another.previous_target_evaluated = 0.0;
  another.target_evaluated = 0.0;
  another.subsample_previous_target_evaluated = 0.0;
  another.subsample_target_evaluated = 0.0;
  another.previous_ensemble_target_evaluated = 0.0;
  another.ensemble_target_evaluated = 0.0;
  another.subsample_previous_ensemble_target_evaluated = 0.0;
  another.subsample_ensemble_target_evaluated = 0.0;
  another.parameters = Parameters();
  another.gradient_estimator_outputs = boost::unordered_map<const ProposalKernel*, GradientEstimatorOutput*>();
  another.move_transformed_parameters = Parameters();
  another.move_transform = NULL;
  another.move_parameters = NULL;
  another.accepted_outputs = boost::unordered_map< const ProposalKernel*, bool>();
  another.previous_self = NULL;
}

Particle Particle::copy_without_factor_variables() const
{
  Particle new_particle;
  
  new_particle.previous_target_evaluated = this->previous_target_evaluated;
  new_particle.target_evaluated = this->target_evaluated;
  new_particle.subsample_previous_target_evaluated = this->subsample_previous_target_evaluated;
  new_particle.subsample_target_evaluated = this->subsample_target_evaluated;
  
  new_particle.previous_ensemble_target_evaluated = this->previous_ensemble_target_evaluated;
  new_particle.subsample_previous_ensemble_target_evaluated = this->subsample_previous_ensemble_target_evaluated;
  new_particle.ensemble_target_evaluated = this->ensemble_target_evaluated;
  new_particle.subsample_ensemble_target_evaluated = this->subsample_ensemble_target_evaluated;
  
  new_particle.parameters = this->parameters;
  new_particle.parameters.self_deep_copy_nonfixed();
  
  new_particle.gradient_estimator_outputs = boost::unordered_map<const ProposalKernel*, GradientEstimatorOutput*>();
  new_particle.gradient_estimator_outputs.reserve(this->gradient_estimator_outputs.size());
  for (auto i=this->gradient_estimator_outputs.begin();
       i!=this->gradient_estimator_outputs.end();
       ++i)
  {
    if (i->second!=NULL)
      new_particle.gradient_estimator_outputs[i->first] = i->second->duplicate();
    else
      new_particle.gradient_estimator_outputs[i->first] = NULL;
  }
  
  new_particle.move_transformed_parameters = this->move_transformed_parameters;
  new_particle.move_transform = this->move_transform;
  
  if (new_particle.move_parameters==&this->parameters)
  {
    new_particle.move_parameters = &new_particle.parameters;
  }
  else if (new_particle.move_parameters==&this->move_transformed_parameters)
  {
    new_particle.move_parameters = &new_particle.move_transformed_parameters;
  }
  
  new_particle.accepted_outputs = this->accepted_outputs;
  
  new_particle.previous_self = this->previous_self;
  
  return new_particle;
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

/*
void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                     const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index,
                                                                  conditioned_on_parameters);
}
*/

void Particle::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->subsample_evaluate_smcfixed_part_of_likelihoods(index);
}

/*
void Particle::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->subsample_evaluate_smcfixed_part_of_likelihoods(index,
                                                                            conditioned_on_parameters);
}
*/

double Particle::evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  //return this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
  
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

/*
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
*/

double Particle::subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->factor_variables->subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    return this->target_evaluated;
  }
  else
    return 0.0;
}

/*
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
*/

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

/*
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
*/

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

/*
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

/*
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
*/

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

/*
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
*/

double Particle::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature);
  else
    return 0.0;
}

/*
double Particle::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature,
                                                     const Parameters &conditioned_on_parameters)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature,
                                                                                conditioned_on_parameters);
  else
    return 0.0;
}
*/

double Particle::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature);
  else
    return 0.0;
  
}

/*
double Particle::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature,
                                                               const Parameters &conditioned_on_parameters)
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature,
                                                                                conditioned_on_parameters);
  else
    return 0.0;
  
}
*/

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
    this->ensemble_target_evaluated = this->ensemble_factor_variables->get_ensemble_factors()->get_temperature()*this->ensemble_target_evaluated;
  }
  return this->target_evaluated + this->ensemble_target_evaluated;
}

/*
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
    this->ensemble_target_evaluated = this->ensemble_factor_variables->get_ensemble_factors()->get_temperature()*this->ensemble_target_evaluated;
  }
  return this->target_evaluated + this->ensemble_target_evaluated;
}
*/

double Particle::subsample_evaluate_all_likelihoods(const Index* index)
{
  this->subsample_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
  if (this->factor_variables!=NULL)
  {
    this->subsample_target_evaluated = this->subsample_target_evaluated + this->factor_variables->subsample_evaluate_likelihoods(index);
  }
  if (this->ensemble_factor_variables!=NULL)
  {
    this->subsample_ensemble_target_evaluated = this->subsample_ensemble_target_evaluated + this->ensemble_factor_variables->subsample_evaluate_likelihoods(index);
  }
  return this->subsample_target_evaluated + this->ensemble_factor_variables->get_ensemble_factors()->get_temperature()*this->subsample_ensemble_target_evaluated;
}

/*
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
  return this->subsample_target_evaluated + this->ensemble_factor_variables->get_ensemble_factors()->get_temperature()*this->subsample_ensemble_target_evaluated;
}
*/

//arma::colvec Particle::get_vector() const
//{
//  return this->parameters.get_vector();
//}

arma::colvec Particle::get_colvec(const std::string &state_name) const
{
  return this->parameters.get_colvec(state_name);
}

arma::rowvec Particle::get_rowvec(const std::string &state_name) const
{
  return this->parameters.get_rowvec(state_name);
}

arma::colvec Particle::get_colvec(const std::vector<std::string> &state_names) const
{
  return this->parameters.get_colvec(state_names);
}

arma::rowvec Particle::get_rowvec(const std::vector<std::string> &state_names) const
{
  return this->parameters.get_rowvec(state_names);
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

void Particle::simulate_factor_variables(Factors* factors)
{
  FactorVariables* new_factor_variables = factors->simulate_factor_variables(this->parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    delete this->factor_variables;
  }
  
  this->factor_variables = new_factor_variables;
}

void Particle::simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors)
{
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->simulate_ensemble_factor_variables(this->parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    delete this->ensemble_factor_variables;
  }
  
  this->ensemble_factor_variables = new_ensemble_factor_variables;
}

/*
void Particle::simulate_factor_variables(Factors* factors)
{
  FactorVariables* new_factor_variables = factors->simulate_factor_variables(this->parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    delete this->factor_variables;
  }
  
  this->factor_variables = new_factor_variables;
}
*/

/*
void Particle::simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors,
                                                  const Parameters &conditioned_on_parameters)
{
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->simulate_ensemble_factor_variables(this->parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    delete this->ensemble_factor_variables;
  }
  
  this->ensemble_factor_variables = new_ensemble_factor_variables;
}
*/

void Particle::subsample_simulate_factor_variables(Factors* factors)
{
  FactorVariables* new_factor_variables = factors->subsample_simulate_factor_variables(this->parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    delete this->factor_variables;
  }
  
  this->factor_variables = new_factor_variables;
}

void Particle::subsample_simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors)
{
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->subsample_simulate_ensemble_factor_variables(this->parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    delete this->ensemble_factor_variables;
  }
  
  this->ensemble_factor_variables = new_ensemble_factor_variables;
}

/*
void Particle::simulate_factor_variables(Factors* factors,
                                         const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  
  FactorVariables* new_factor_variables = factors->simulate_factor_variables(all_parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    delete this->factor_variables;
  }
  
  this->factor_variables = new_factor_variables;
}

void Particle::simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors,
                                                  const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->simulate_ensemble_factor_variables(all_parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    delete this->ensemble_factor_variables;
  }
  
  this->ensemble_factor_variables = new_ensemble_factor_variables;
}

void Particle::subsample_simulate_factor_variables(Factors* factors,
                                                   const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  
  FactorVariables* new_factor_variables = factors->subsample_simulate_factor_variables(all_parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    delete this->factor_variables;
  }
  
  this->factor_variables = new_factor_variables;
}

void Particle::subsample_simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors,
                                                            const Parameters &conditioned_on_parameters)
{
  Parameters all_parameters = this->parameters.merge(conditioned_on_parameters);
  
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->subsample_simulate_ensemble_factor_variables(all_parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    delete this->ensemble_factor_variables;
  }
  
  this->ensemble_factor_variables = new_ensemble_factor_variables;
}
*/

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
