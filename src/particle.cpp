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

namespace ilike
{
Particle::Particle()
: parameters()
//move_transformed_parameters()
{
  
  this->factor_variables.reset();
  this->ensemble_factor_variables.reset();
  
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

Particle::Particle(const Parameters &parameters_in)
: parameters(parameters_in)
{
  
  this->factor_variables.reset();
  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
: parameters(parameters_in)
//move_transformed_parameters()
{
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                   const Parameters &conditioned_on_parameters)
: parameters(parameters_in)
//move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
: parameters(parameters_in)
//move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
                   FactorVariables* factor_variables_in)
: parameters(parameters_in),
//move_transformed_parameters(),
factor_variables(factor_variables_in)
{
  
  this->ensemble_factor_variables.reset();
  
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

Particle::Particle(const Parameters &parameters_in,
                   FactorVariables* factor_variables_in,
                   double previous_target_evaluated_in)
: parameters(parameters_in),
factor_variables(factor_variables_in)
{
  
  this->ensemble_factor_variables.reset();
  
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

Particle::Particle(const Parameters &parameters_in,
                   const EnsembleFactors* ensemble_factors_in)
: parameters(parameters_in)
//move_transformed_parameters()
{
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }

  this->factor_variables.reset();
  
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
                   const EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters)
{
  this->parameters = parameters_in;
  
  this->parameters.merge_with_fixed(conditioned_on_parameters);

  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables.reset();
  }
  this->factor_variables.reset();
  
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
                   const EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
{
  this->parameters = parameters_in;
  
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);

  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables.reset();
  }
  this->factor_variables.reset();
  
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
                   EnsembleFactorVariables* ensemble_factor_variables_in)
{
  this->parameters = parameters_in;
  
  this->ensemble_factor_variables.reset(ensemble_factor_variables_in);
  this->factor_variables.reset();
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

Particle::Particle(Parameters &&parameters_in)
: parameters(std::move(parameters_in))
{
  
  this->factor_variables.reset();
  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
: parameters(std::move(parameters_in))
//move_transformed_parameters()
{
  
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                   const Parameters &conditioned_on_parameters)
: parameters(std::move(parameters_in))
//move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
                   const Factors* factors_in,
                   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
: parameters(std::move(parameters_in))
//move_transformed_parameters()
{
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);
  
  if (factors_in!=NULL)
  {
    this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
  
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);

  this->ensemble_factor_variables.reset();
  
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
//move_transformed_parameters(),
factor_variables(factor_variables_in)
{
  
  this->ensemble_factor_variables.reset();
  
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
factor_variables(factor_variables_in)
{
  
  this->ensemble_factor_variables.reset();
  
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
                   const EnsembleFactors* ensemble_factors_in)
: parameters(std::move(parameters_in))
//move_transformed_parameters()
{
  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }

  this->factor_variables.reset();
  
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
                   const EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters)
{
  this->parameters = std::move(parameters_in);
  
  this->parameters.merge_with_fixed(conditioned_on_parameters);

  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables.reset();
  }
  this->factor_variables.reset();
  
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
                   const EnsembleFactors* ensemble_factors_in,
                   const Parameters &conditioned_on_parameters,
                   const Parameters &sequencer_parameters)
{
  this->parameters = std::move(parameters_in);
  
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);

  if (ensemble_factors_in!=NULL)
  {
    this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }
  else
  {
    this->ensemble_factor_variables.reset();
  }
  this->factor_variables.reset();
  
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
  
  this->ensemble_factor_variables.reset(ensemble_factor_variables_in);
  this->factor_variables.reset();
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
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(Factors* factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(Factors* factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);
  
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in)
{
  this->parameters = parameters_in;
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     Factors* factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);
  
  this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  this->factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in)
{
  this->parameters = parameters_in;
  this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
  this->ensemble_factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in,
                     const Parameters &conditioned_on_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters);
  
  this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
  this->ensemble_factor_variables->set_particle(this);
}

void Particle::setup(const Parameters &parameters_in,
                     EnsembleFactors* ensemble_factors_in,
                     const Parameters &conditioned_on_parameters,
                     const Parameters &sequencer_parameters)
{
  this->parameters = parameters_in;
  this->parameters.merge_with_fixed(conditioned_on_parameters,
                                    sequencer_parameters);
  
  this->ensemble_factor_variables.reset(ensemble_factors_in->simulate_ensemble_factor_variables(this->parameters));
  this->ensemble_factor_variables->set_particle(this);
}

Particle::~Particle()
{
  this->factor_variables.reset();
  
  this->ensemble_factor_variables.reset();
  
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
  
  this->factor_variables.reset();
  
  this->ensemble_factor_variables.reset();

  this->make_copy(another);
  
  return *this;
}

void Particle::make_copy(const Particle &another)
{
  if (another.factor_variables!=NULL)
  {
    this->factor_variables.reset(another.factor_variables->duplicate());
    this->factor_variables->set_particle(this);
  }
  else
    this->factor_variables.reset();
  
  if (another.ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables.reset(another.ensemble_factor_variables->duplicate());
    this->ensemble_factor_variables->set_particle(this);
  }
  else
    this->ensemble_factor_variables.reset();
  
  this->previous_target_evaluated = another.previous_target_evaluated;
  this->target_evaluated = another.target_evaluated;
  this->subsample_previous_target_evaluated = another.subsample_previous_target_evaluated;
  this->subsample_target_evaluated = another.subsample_target_evaluated;
  
  this->previous_ensemble_target_evaluated = another.previous_ensemble_target_evaluated;
  this->subsample_previous_ensemble_target_evaluated = another.subsample_previous_ensemble_target_evaluated;
  this->ensemble_target_evaluated = another.ensemble_target_evaluated;
  this->subsample_ensemble_target_evaluated = another.subsample_ensemble_target_evaluated;
  this->proposal_serial = another.proposal_serial;
  this->proposal_log_ratio = another.proposal_log_ratio;
  
  this->parameters = another.parameters;

  this->accepted_outputs = another.accepted_outputs;
  
  this->previous_self = another.previous_self;
  
  this->current_proposal_store = another.current_proposal_store;
  
  this->proposals_to_transform_for_pointer = another.proposals_to_transform_for_pointer;
  this->proposals_to_find_gradient_for_pointer = another.proposals_to_find_gradient_for_pointer;
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
  
  this->factor_variables.reset();
  
  this->ensemble_factor_variables.reset();

  this->make_copy(std::move(another));
  
  return *this;
}

void Particle::make_copy(Particle &&another)
{
  if (another.factor_variables!=NULL)
  {
    this->factor_variables = std::move(another.factor_variables);
    this->factor_variables->set_particle(this);
  }
  else
    this->factor_variables.reset();
  
  if (another.ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables = std::move(another.ensemble_factor_variables);
    this->ensemble_factor_variables->set_particle(this);
  }
  else
    this->ensemble_factor_variables.reset();
  
  this->previous_target_evaluated = std::move(another.previous_target_evaluated);
  this->target_evaluated = std::move(another.target_evaluated);
  this->proposal_serial = another.proposal_serial;
  this->proposal_log_ratio = another.proposal_log_ratio;
  this->subsample_previous_target_evaluated = std::move(another.subsample_previous_target_evaluated);
  this->subsample_target_evaluated = std::move(another.subsample_target_evaluated);
  
  this->previous_ensemble_target_evaluated = std::move(another.previous_ensemble_target_evaluated);
  this->subsample_previous_ensemble_target_evaluated = std::move(another.subsample_previous_ensemble_target_evaluated);
  this->ensemble_target_evaluated = std::move(another.ensemble_target_evaluated);
  this->subsample_ensemble_target_evaluated = std::move(another.subsample_ensemble_target_evaluated);
  
  this->parameters = std::move(another.parameters);

  this->accepted_outputs = std::move(another.accepted_outputs);
  
  this->previous_self = another.previous_self;
  
  //this->current_proposal_store = boost::unordered_map<const ProposalKernel*, ProposalStore>();
  // Move the whole map (its nodes change owner; no element is copied).
  this->current_proposal_store = std::move(another.current_proposal_store);
  
  this->proposals_to_transform_for_pointer = another.proposals_to_transform_for_pointer;
  this->proposals_to_find_gradient_for_pointer = another.proposals_to_find_gradient_for_pointer;

  another.factor_variables.reset();
  another.ensemble_factor_variables.reset();
  another.previous_target_evaluated = 0.0;
  another.target_evaluated = 0.0;
  another.subsample_previous_target_evaluated = 0.0;
  another.subsample_target_evaluated = 0.0;
  another.previous_ensemble_target_evaluated = 0.0;
  another.ensemble_target_evaluated = 0.0;
  another.subsample_previous_ensemble_target_evaluated = 0.0;
  another.subsample_ensemble_target_evaluated = 0.0;
  // The moved-from containers are emptied with clear(), not by assigning
  // newly constructed ones (a temporary built and destroyed per move).
  // another.parameters was emptied by its own move above.
  another.accepted_outputs.clear();
  another.previous_self = NULL;
  another.current_proposal_store.clear();
  
  another.proposals_to_transform_for_pointer = NULL;
  another.proposals_to_find_gradient_for_pointer = NULL;
}

void Particle::simulate_factor_variables()
{
  if (this->factor_variables->get_factors()!=NULL)
  {
    this->factor_variables.reset(this->factor_variables->get_factors()->simulate_factor_variables(this->parameters));
    this->factor_variables->set_particle(this);
  }
}

void Particle::simulate_ensemble_factor_variables()
{
  if (this->ensemble_factor_variables->get_ensemble_factors()!=NULL)
  {
    this->ensemble_factor_variables.reset(this->ensemble_factor_variables->get_ensemble_factors()->simulate_ensemble_factor_variables(this->parameters));
    this->ensemble_factor_variables->set_particle(this);
  }
}

void Particle::refill(Parameters &&parameters_in,
                      const Factors* factors_in,
                      const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                      const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
                      const Parameters &conditioned_on_parameters,
                      const Parameters* second_conditioned_on_parameters)
{
  // as the constructors: merge_with_fixed(conditioned) or
  // merge_with_fixed(conditioned, sequencer)
  this->parameters = std::move(parameters_in);
  if (second_conditioned_on_parameters==NULL)
    this->parameters.merge_with_fixed(conditioned_on_parameters);
  else
    this->parameters.merge_with_fixed(conditioned_on_parameters, *second_conditioned_on_parameters);

  bool recycled = false;
  if ((factors_in!=NULL) && (this->factor_variables!=NULL) &&
      (this->factor_variables->get_factors()==factors_in))
    recycled = this->factor_variables->recycle_all_outputs(factors_in, this->parameters);
  if (!recycled)
  {
    this->factor_variables.reset();
    this->factor_variables.reset();
    if (factors_in!=NULL)
      this->factor_variables.reset(factors_in->simulate_factor_variables(this->parameters));
  }
  if (this->factor_variables!=NULL)
    this->factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables.reset();
    this->ensemble_factor_variables.reset();
  }

  // as a new Particle
  this->current_proposal_store.clear();
  this->accepted_outputs.clear();
  this->previous_self = NULL;
  this->simulate_proposal_variables(proposals_to_transform_for_in, proposals_to_find_gradient_for_in);
  this->previous_target_evaluated = 0.0;
  this->target_evaluated = 0.0;
  this->subsample_previous_target_evaluated = 0.0;
  this->subsample_target_evaluated = 0.0;
  this->previous_ensemble_target_evaluated = 0.0;
  this->subsample_previous_ensemble_target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  this->subsample_ensemble_target_evaluated = 0.0;
}

void Particle::simulate_proposal_variables(const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                                           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in)
{
  this->proposals_to_find_gradient_for_pointer = proposals_to_find_gradient_for_in;
  this->proposals_to_transform_for_pointer = proposals_to_transform_for_in;
  
  if (this->proposals_to_find_gradient_for_pointer!=NULL)
  {
    for (auto i=this->proposals_to_find_gradient_for_pointer->begin();
         i!=this->proposals_to_find_gradient_for_pointer->end();
         ++i)
    {
      GradientEstimatorOutput* simulated = (*i)->simulate_gradient_estimator_output();
      if (simulated!=NULL)
      {
        this->current_proposal_store[(*i)->get_instance_index()] = ProposalStore(simulated);
      }
    }
  }
  
  if (this->proposals_to_transform_for_pointer!=NULL)
  {
    for (auto i=this->proposals_to_transform_for_pointer->begin();
         i!=proposals_to_transform_for_pointer->end();
         ++i)
    {
      // A proposal without a transform needs nothing stored:
      // get_transformed_parameters returns this->parameters for it without
      // looking in the store. Storing a copy of the parameters here cost a
      // hash-map copy (and a store entry) per proposed particle.
      if ((*i)->get_transform()==NULL)
        continue;

      boost::unordered_map< int, ProposalStore>::iterator current_transform_stored_at_proposal = this->current_proposal_store.end();

      auto found = this->current_proposal_store.find((*i)->get_instance_index());
      
      if (found != this->current_proposal_store.end()) // if proposal is found in current_proposal_store
      {
        
        if ((*i)->get_transform()==NULL)
        {
          found->second.set_transformed_parameters(this->parameters);
        }
        else
        {
          found->second.set_transformed_parameters((*i)->get_transform()->transform(this->parameters));
        }
        
      }
      else
      {
        
        if ((*i)->get_transform()==NULL)
        {
          this->current_proposal_store[(*i)->get_instance_index()] = this->parameters;
        }
        else
        {
          this->current_proposal_store[(*i)->get_instance_index()] = (*i)->get_transform()->transform(this->parameters);
        }
        
      }
    }
    
  }
  
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

  new_particle.accepted_outputs = this->accepted_outputs;
  
  new_particle.previous_self = this->previous_self;
  
  return new_particle;
}

void Particle::evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
  {
    this->factor_variables->evaluate_smcfixed_part_of_likelihoods(index);
  }
}

void Particle::subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)
{
  if (this->factor_variables!=NULL)
    this->factor_variables->subsample_evaluate_smcfixed_part_of_likelihoods(index);
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

double Particle::get_cached_likelihoods(const Index* index) const
{
  if (this->factor_variables!=NULL)
    return this->factor_variables->get_cached_likelihoods(index);
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

arma::mat Particle::direct_get_gradient_of_log(const std::string &variable,
                                               const Index* index) const
{
  // a stacked variable (C203): the stack of its variables' gradients
  // (C204: each factor is asked once for all of its variables)
  if (is_stacked_name(variable))
  {
    std::vector<std::string> names = split_stacked_name(variable);
    if (this->factor_variables==NULL)
      return arma::zeros<arma::colvec>(this->parameters.get_colvec(variable).n_elem);
    std::vector<arma::mat> gradients = this->factor_variables->direct_get_gradients_of_log(index, names);
    arma::uword total = 0;
    for (const arma::mat &gradient : gradients)
      total += gradient.n_elem;
    arma::colvec stacked(total);
    arma::uword offset = 0;
    for (const arma::mat &gradient : gradients)
    {
      if (gradient.n_elem > 0)
        stacked.subvec(offset, offset + gradient.n_elem - 1) = arma::vectorise(gradient);
      offset += gradient.n_elem;
    }
    return stacked;
  }
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_get_gradient_of_log(index,variable);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

arma::mat Particle::direct_subsample_get_gradient_of_log(const std::string &variable,
                                                         const Index* index) const
{
  // a stacked variable (C203): the stack of its variables' gradients
  if (is_stacked_name(variable))
  {
    arma::colvec stacked;
    for (const std::string &name : split_stacked_name(variable))
    {
      arma::colvec gradient = arma::vectorise(this->direct_subsample_get_gradient_of_log(name, index));
      stacked = (stacked.n_elem == 0) ? gradient : arma::colvec(arma::join_cols(stacked, gradient));
    }
    return stacked;
  }
  arma::mat current_parameter = this->parameters[variable];
  
  if (this->factor_variables!=NULL)
  {
    return this->factor_variables->direct_subsample_get_gradient_of_log(index,
                                                                        variable);
  }
  else
    return arma::mat(current_parameter.n_rows,current_parameter.n_cols);
}

double Particle::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature,
                                                     const std::vector<arma::mat> &inv_sigma_precomps,
                                                     const std::vector<double> &log_det_precomps) const
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature,
                                                                                inv_sigma_precomps,
                                                                                log_det_precomps);
  else
    return 0.0;
}

double Particle::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature,
                                                               const std::vector<arma::mat> &inv_sigma_precomps,
                                                               const std::vector<double> &log_det_precomps) const
{
  if (this->ensemble_factor_variables!=NULL)
    return this->ensemble_factor_variables->evaluate_ensemble_likelihood_ratios(index,
                                                                                inverse_incremental_temperature,
                                                                                inv_sigma_precomps,
                                                                                log_det_precomps);
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
    
    Rcpp::stop("Particle::evaluate_all_likelihoods - need temperature to calculate ensemble llhds.");
  }
  return this->target_evaluated + this->ensemble_target_evaluated;
}

double Particle::subsample_evaluate_all_likelihoods(const Index* index)
{
  this->target_evaluated = 0.0;
  this->ensemble_target_evaluated = 0.0;
  if (this->factor_variables!=NULL)
  {
    this->target_evaluated = this->target_evaluated + this->factor_variables->subsample_evaluate_likelihoods(index);
  }
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_target_evaluated = this->ensemble_target_evaluated + this->ensemble_factor_variables->subsample_evaluate_likelihoods(index);
    Rcpp::stop("Particle::subsample_evaluate_all_likelihoods - need temperature to calculate ensemble llhds.");
  }
  return this->target_evaluated + this->ensemble_target_evaluated;
}

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

void Particle::swap(Particle &another)
{
  if (this == &another)
    return;
  std::swap(this->parameters, another.parameters);
  std::swap(this->current_proposal_store, another.current_proposal_store);
  std::swap(this->proposals_to_transform_for_pointer, another.proposals_to_transform_for_pointer);
  std::swap(this->proposals_to_find_gradient_for_pointer, another.proposals_to_find_gradient_for_pointer);
  std::swap(this->accepted_outputs, another.accepted_outputs);
  std::swap(this->previous_target_evaluated, another.previous_target_evaluated);
  std::swap(this->subsample_previous_target_evaluated, another.subsample_previous_target_evaluated);
  std::swap(this->target_evaluated, another.target_evaluated);
  std::swap(this->subsample_target_evaluated, another.subsample_target_evaluated);
  std::swap(this->previous_ensemble_target_evaluated, another.previous_ensemble_target_evaluated);
  std::swap(this->subsample_previous_ensemble_target_evaluated, another.subsample_previous_ensemble_target_evaluated);
  std::swap(this->ensemble_target_evaluated, another.ensemble_target_evaluated);
  std::swap(this->subsample_ensemble_target_evaluated, another.subsample_ensemble_target_evaluated);
  std::swap(this->proposal_serial, another.proposal_serial);
  std::swap(this->proposal_log_ratio, another.proposal_log_ratio);
  std::swap(this->factor_variables, another.factor_variables);
  std::swap(this->ensemble_factor_variables, another.ensemble_factor_variables);
  std::swap(this->previous_self, another.previous_self);
  if (this->factor_variables!=NULL)
    this->factor_variables->set_particle(this);
  if (another.factor_variables!=NULL)
    another.factor_variables->set_particle(&another);
  if (this->ensemble_factor_variables!=NULL)
    this->ensemble_factor_variables->set_particle(this);
  if (another.ensemble_factor_variables!=NULL)
    another.ensemble_factor_variables->set_particle(&another);
}

const Parameters& Particle::get_transformed_parameters(const ProposalKernel* proposal_in) const
{
  if (proposal_in->get_transform()==NULL)
    return this->parameters;

  auto found = this->current_proposal_store.find(proposal_in->get_instance_index());
  if (found != this->current_proposal_store.end())
  {
    const Parameters &transformed_parameters = found->second.get_transformed_parameters();
    if (!transformed_parameters.is_empty())
      return transformed_parameters;
    found->second.set_transformed_parameters(proposal_in->get_transform()->transform(this->parameters));
    return found->second.get_transformed_parameters();
  }
  ProposalStore &store = this->current_proposal_store[proposal_in->get_instance_index()];
  store.set_transformed_parameters(proposal_in->get_transform()->transform(this->parameters));
  return store.get_transformed_parameters();
}

GradientEstimatorOutput* Particle::get_gradient_estimator_output(const ProposalKernel* proposal_in) const
{
  auto found = this->current_proposal_store.find(proposal_in->get_instance_index());
  if (found != this->current_proposal_store.end())
  {
    return found->second.get_gradient_estimator_output();
  }
  else
  {
    Rcpp::stop("Particle::get_gradient_estimator_output - gradient information not found for proposal.");
  }
}

void Particle::set_acceptance(const ProposalKernel* proposal_in,
                              bool accepted_in)
{
  this->accepted_outputs[proposal_in->get_instance_index()] = accepted_in;
}

void Particle::erase_mcmc_adaptation_info()
{
  this->accepted_outputs.clear();
}

void Particle::simulate_factor_variables(const FactorVariables* current,
                                         const Index* new_outputs_index)
{
  FactorVariables* new_factor_variables = current->copy_with_new_outputs(this->parameters,
                                                                         new_outputs_index);
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    this->factor_variables.reset();
  }
  this->factor_variables.reset(new_factor_variables);
}

void Particle::simulate_factor_variables(const Factors* factors)
{
  FactorVariables* new_factor_variables = factors->simulate_factor_variables(this->parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    this->factor_variables.reset();
  }
  
  this->factor_variables.reset(new_factor_variables);
}

void Particle::simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors)
{
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->simulate_ensemble_factor_variables(this->parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables.reset();
  }
  
  this->ensemble_factor_variables.reset(new_ensemble_factor_variables);
}

void Particle::subsample_simulate_factor_variables(const Factors* factors)
{
  FactorVariables* new_factor_variables = factors->subsample_simulate_factor_variables(this->parameters);
  
  new_factor_variables->set_particle(this);
  if (this->factor_variables!=NULL)
  {
    this->factor_variables.reset();
  }
  
  this->factor_variables.reset(new_factor_variables);
}

void Particle::subsample_simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors)
{
  EnsembleFactorVariables* new_ensemble_factor_variables = ensemble_factors->subsample_simulate_ensemble_factor_variables(this->parameters);
  
  new_ensemble_factor_variables->set_particle(this);
  if (this->ensemble_factor_variables!=NULL)
  {
    this->ensemble_factor_variables.reset();
  }
  
  this->ensemble_factor_variables.reset(new_ensemble_factor_variables);
}

void Particle::tell_factors_to_forget_they_were_already_written_to_file()
{
  if (this->factor_variables!=NULL)
    this->factor_variables->forget_you_were_already_written_to_file();
  
  if (this->ensemble_factor_variables!=NULL)
    this->ensemble_factor_variables->forget_you_were_already_written_to_file();
}

std::ostream& operator<<(std::ostream& os, const Particle &p)
{
  os << p.parameters << std::endl;
  
  os << "{" << std::endl;

  os << "}";
  
  return os;
}
}
