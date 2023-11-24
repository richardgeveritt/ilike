#include <math.h>
#include "proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "composite_proposal_kernel.h"
#include "transform.h"
#include "mcmc_adaptor.h"
#include "smc_adaptor.h"
#include "direct_gradient_estimator_output.h"
#include "factor_variables.h"
#include "factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"

int ProposalKernel::instance_counter = 0;

ProposalKernel::ProposalKernel()
  :Kernel()
{
  this->instance_index = this->instance_counter++;
  this->mcmc_adaptor = NULL;
  this->smc_adaptor = NULL;
  this->transform = NULL;
}

ProposalKernel::~ProposalKernel()
{
  if (this->mcmc_adaptor!=NULL)
    delete this->mcmc_adaptor;
  
  if (this->smc_adaptor!=NULL)
    delete this->smc_adaptor;
  
  if (this->transform!=NULL)
    delete this->transform;
}

ProposalKernel::ProposalKernel(const Parameters &proposal_parameters_in)
  :Kernel()
{
  this->instance_index = this->instance_counter++;
  //this->proposal_parameters = proposal_parameters_in;
}

ProposalKernel::ProposalKernel(const ProposalKernel &another)
  :Kernel(another)
{
  this->make_copy(another);
}

void ProposalKernel::operator=(const ProposalKernel &another)
{
  if(this == &another)
    return;
  
  if (this->mcmc_adaptor!=NULL)
    delete this->mcmc_adaptor;
  
  if (this->smc_adaptor!=NULL)
    delete this->smc_adaptor;
  
  if (this->transform!=NULL)
    delete this->transform;

  Kernel::operator=(another);
  this->make_copy(another);
}

void ProposalKernel::make_copy(const ProposalKernel &another)
{
  if (another.mcmc_adaptor!=NULL)
    this->mcmc_adaptor = another.mcmc_adaptor->duplicate();
  else
    this->mcmc_adaptor = NULL;
  
  if (another.smc_adaptor!=NULL)
    this->smc_adaptor = another.smc_adaptor->duplicate();
  else
    this->smc_adaptor = NULL;
  
  if (another.transform!=NULL)
    this->transform = another.transform->duplicate();
  else
    this->transform = NULL;
  
  //this->proposal_parameters = another.proposal_parameters;
  
  this->instance_index = another.instance_index;
  this->instance_counter = another.instance_counter;
}

Particle ProposalKernel::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  Particle proposed_particle;
  proposed_particle.parameters = particle.parameters;
  
  if (this->transform==NULL)
  {
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->simulate(rng,particle));
  }
  else
  {
    Parameters proposed_parameters_in_transformed_space = this->simulate(rng,particle);
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  
  // Outputs are created here, with memory managed by Particle hereafter.
  if (particle.factor_variables!=NULL)
  {
    const Factors* old_factors = particle.factor_variables->get_factors();
    if (old_factors!=NULL)
      proposed_particle.simulate_factor_variables(old_factors);
  }
  
  if (particle.ensemble_factor_variables!=NULL)
  {
    const EnsembleFactors* old_ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (old_ensemble_factors!=NULL)
      proposed_particle.simulate_ensemble_factor_variables(old_ensemble_factors);
  }
  
  proposed_particle.simulate_proposal_variables(particle.proposals_to_transform_for_pointer, particle.proposals_to_find_gradient_for_pointer);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  proposed_particle.previous_self = &particle;
  
  return proposed_particle; // what we need is for previous target eval to be set to eval for this target (when evaluated there)!!!!! we don't want it to know its grad, until evaluated there
}

/*
Particle ProposalKernel::move(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->simulate(rng,
                                                  particle,
                                                  conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->simulate(rng,
                                                                   particle,
                                                                   conditioned_on_parameters);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  if (particle.factor_variables!=NULL)
  {
    Factors* old_factors = particle.factor_variables->get_factors();
    if (old_factors!=NULL)
      proposed_particle.simulate_factor_variables(old_factors,
                                                  conditioned_on_parameters);
  }
  
  if (particle.ensemble_factor_variables!=NULL)
  {
    EnsembleFactors* old_ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (old_ensemble_factors!=NULL)
      proposed_particle.simulate_ensemble_factor_variables(old_ensemble_factors,
                                                           conditioned_on_parameters);
  }
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  proposed_particle.previous_self = &particle;
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}
*/

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  Particle proposed_particle;
  proposed_particle.parameters = particle.parameters;
  
  if (this->transform==NULL)
  {
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->subsample_simulate(rng,particle));
  }
  else
  {
    Parameters proposed_parameters_in_transformed_space = this->subsample_simulate(rng,particle);
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  
  // Outputs are created here, with memory managed by Particle hereafter.
  if (particle.factor_variables!=NULL)
  {
    const Factors* old_factors = particle.factor_variables->get_factors();
    if (old_factors!=NULL)
      proposed_particle.subsample_simulate_factor_variables(old_factors);
  }
  
  if (particle.ensemble_factor_variables!=NULL)
  {
    const EnsembleFactors* old_ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (old_ensemble_factors!=NULL)
      proposed_particle.subsample_simulate_ensemble_factor_variables(old_ensemble_factors);
  }
  
  proposed_particle.simulate_proposal_variables(particle.proposals_to_transform_for_pointer, particle.proposals_to_find_gradient_for_pointer);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  proposed_particle.previous_self = &particle;
  
  return proposed_particle; // what we need is for previous target eval to be set to eval for this target (when evaluated there)!!!!! we don't want it to know its grad, until evaluated there
}

/*
Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,particle,conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,particle);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  proposed_particle.factor_variables = particle.factor_variables->get_factors()->simulate_factor_variables(all_proposed_parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  
  proposed_particle.previous_self = &particle;
  
  return proposed_particle;
}
*/

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        const Particle &particle) const
{
  Particle proposed_particle;
  proposed_particle.parameters = particle.parameters;
  
  if (this->transform==NULL)
  {
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->subsample_simulate(rng,variable,particle));
  }
  else
  {
    Parameters proposed_parameters_in_transformed_space = this->subsample_simulate(rng,variable,particle);
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  
  // Outputs are created here, with memory managed by Particle hereafter.
  if (particle.factor_variables!=NULL)
  {
    const Factors* old_factors = particle.factor_variables->get_factors();
    if (old_factors!=NULL)
      proposed_particle.subsample_simulate_factor_variables(old_factors);
  }
  
  if (particle.ensemble_factor_variables!=NULL)
  {
    const EnsembleFactors* old_ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (old_ensemble_factors!=NULL)
      proposed_particle.subsample_simulate_ensemble_factor_variables(old_ensemble_factors);
  }
  
  proposed_particle.simulate_proposal_variables(particle.proposals_to_transform_for_pointer, particle.proposals_to_find_gradient_for_pointer);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  //proposed_particle.previous_self = &particle;
  
  return proposed_particle; // what we need is for previous target eval to be set to eval for this target (when evaluated there)!!!!! we don't want it to know its grad, until evaluated there
}

/*
Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,
                                                            variable,
                                                            particle,
                                                            conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,
                                                                             variable,
                                                                             particle,
                                                                             conditioned_on_parameters);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  proposed_particle.factor_variables = particle.factor_variables->get_factors()->subsample_simulate_factor_variables(all_proposed_parameters,
                                                                                                                     conditioned_on_parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  proposed_particle.previous_self = &particle;
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}
*/

/*
Particle ProposalKernel::move(RandomNumberGenerator &rng,
                                    const Index* index,
                                    Particle &particle) const
{
  
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->simulate(rng,particle);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->simulate(rng,particle);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
    
    // thing to think about...
    // reason we pass particle is because it has extra info such as gradient
    // we sometimes set this in the proposal (such as finding grad)
    // if we pass temp particle, this info won't be stored...
    // what info do we want stored? gradient of transformed? not sure...
    // might need to store transformed version in Particle? Try to avoid
    
    // put in simulate!
    // deriv in some algs comes from deriv of prior times llhd in transformed space - need diff bit
  }
  
  // Outputs are created here, with memory managed by Particle hereafter.
  //if (particle.ensemble_factor_variables!=NULL)
  //  proposed_particle.ensemble_factor_variables = particle.ensemble_factor_variables->ensemble_factors->simulate_ensemble_factor_variables(rng,
                                                                                                                                         proposed_particle.parameters);
  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->get_factors()->simulate_factor_variables(rng,
                                                                                                             proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  
  return proposed_particle; // what we need is for previous target eval to be set to eval for this target (when evaluated there)!!!!! we don't want it to know its grad, until evaluated there
}

Particle ProposalKernel::move(RandomNumberGenerator &rng,
                                    const Index* index,
                                    Particle &particle,
                                    const Parameters &conditioned_on_parameters) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->simulate(rng,particle,conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->simulate(rng,
                                                                   particle,
                                                                   conditioned_on_parameters);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  //Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  //if (particle.ensemble_factor_variables!=NULL)
  //  proposed_particle.ensemble_factor_variables = particle.ensemble_factor_variables->ensemble_factors->simulate_ensemble_factor_variables(rng,
                                                                                                                                         all_proposed_parameters);
  
  
  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->get_factors()->simulate_factor_variables(rng,
                                                                                                             proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                              Particle &particle) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,particle);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,particle);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->factors->simulate_factor_variables(rng,
                                                                                                       proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                              const Index* index,
                                              Particle &particle,
                                              const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,
                                                            particle,
                                                            conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,
                                                                             particle,
                                                                             conditioned_on_parameters);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  //Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  //// Outputs are created here, with memory managed by Particle hereafter.
  //if (particle.ensemble_factor_variables!=NULL)
  //  proposed_particle.ensemble_factor_variables = particle.ensemble_factor_variables->ensemble_factors->simulate_ensemble_factor_variables(rng,
                                                                                                                                         all_proposed_parameters);
  

  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->get_factors()->simulate_factor_variables(rng,
                                                                                                             proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                              const std::string &variable,
                                              Particle &particle) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,variable,particle);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,variable,particle);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->factors->simulate_factor_variables(rng,
                                                                                                       proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}

Particle ProposalKernel::subsample_move(RandomNumberGenerator &rng,
                                              const Index* index,
                                              const std::string &variable,
                                              Particle &particle,
                                              const Parameters &conditioned_on_parameters) const
{
  //Parameters all_parameters = particle.parameters.merge(conditioned_on_parameters);
  Particle proposed_particle;
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
    proposed_particle.parameters = this->subsample_simulate(rng,
                                                            variable,
                                                            particle,
                                                            conditioned_on_parameters);
  }
  else
  {
    // transform particles
    particle.set_move_transformed_parameters(this->transform);
    proposed_particle.move_transformed_parameters = this->subsample_simulate(rng,
                                                                             variable,
                                                                             particle,
                                                                             conditioned_on_parameters);
    proposed_particle.parameters = this->transform->inverse_transform(proposed_particle.move_transformed_parameters);
  }
  
  //Parameters all_proposed_parameters = proposed_particle.parameters.merge(conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  //if (particle.ensemble_factor_variables!=NULL)
  //  proposed_particle.ensemble_factor_variables = particle.ensemble_factor_variables->ensemble_factors->aubsample_simulate_ensemble_factor_variables(rng,
                                                                                                                                                   all_proposed_parameters);
  
  if (particle.factor_variables!=NULL)
    proposed_particle.factor_variables = particle.factor_variables->get_factors()->simulate_factor_variables(rng,
                                                                                                             proposed_particle.parameters);
  
  proposed_particle.accepted_outputs = particle.accepted_outputs;
  // Outputs are created here, with memory managed by Particle hereafter.
  return proposed_particle;
}
*/

void ProposalKernel::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  if (this->smc_adaptor!=NULL)
    this->smc_adaptor->ensemble_adapt(current_state);
}

void ProposalKernel::smc_adapt(SMCOutput* current_state)
{
  if (this->smc_adaptor!=NULL)
    this->smc_adaptor->smc_adapt(current_state);
}

void ProposalKernel::mcmc_adapt(const Particle &current_particle,
                                size_t iteration_counter)
{
  if (this->mcmc_adaptor!=NULL)
    this->mcmc_adaptor->mcmc_adapt(current_particle,
                                   iteration_counter);
}

/*
void ProposalKernel::use_transform(Particle &particle)
{
  if (this->transform==NULL)
  {
    particle.set_move_transformed_parameters();
  }
  else
  {
    particle.set_move_transformed_parameters(this->transform);
  }
}
*/

double ProposalKernel::evaluate_kernel(const Particle &proposed_particle,
                                       const Particle &old_particle) const
{
  if (this->transform==NULL)
  {
    //proposed_particle.set_move_transformed_parameters();
    //old_particle.set_move_transformed_parameters();
    return this->specific_evaluate_kernel(proposed_particle, old_particle);
  }
  else
  {
    //proposed_particle.set_move_transformed_parameters(this->transform);
    //old_particle.set_move_transformed_parameters(this->transform);
    
    // need to include absolute value of Jacobian determinant
    return this->specific_evaluate_kernel(proposed_particle, old_particle) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}

/*
double ProposalKernel::evaluate_kernel(Particle &proposed_particle,
                                       Particle &old_particle,
                                       const Parameters &conditioned_on_parameters) const
{
  if (this->transform==NULL)
  {
    proposed_particle.set_move_transformed_parameters();
    old_particle.set_move_transformed_parameters();
    return this->specific_evaluate_kernel(proposed_particle,
                                          old_particle,
                                          conditioned_on_parameters);
  }
  else
  {
    proposed_particle.set_move_transformed_parameters(this->transform);
    old_particle.set_move_transformed_parameters(this->transform);
    
    // need to include absolute value of Jacobian determinant
    return this->specific_evaluate_kernel(proposed_particle,
                                          old_particle,
                                          conditioned_on_parameters) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}
*/

double ProposalKernel::subsample_evaluate_kernel(const Particle &proposed_particle,
                                                 const Particle &old_particle) const
{
  if (this->transform==NULL)
  {
    //proposed_particle.set_move_transformed_parameters();
    //old_particle.set_move_transformed_parameters();
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle);
  }
  else
  {
    //proposed_particle.set_move_transformed_parameters(this->transform);
    //old_particle.set_move_transformed_parameters(this->transform);
    
    // need to include absolute value of Jacobian determinant
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}

/*
double ProposalKernel::subsample_evaluate_kernel(Particle &proposed_particle,
                                                 Particle &old_particle,
                                                 const Parameters &conditioned_on_parameters) const
{
  if (this->transform==NULL)
  {
    proposed_particle.set_move_transformed_parameters();
    old_particle.set_move_transformed_parameters();
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle,
                                                    conditioned_on_parameters);
  }
  else
  {
    proposed_particle.set_move_transformed_parameters(this->transform);
    old_particle.set_move_transformed_parameters(this->transform);
    
    // need to include absolute value of Jacobian determinant
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle,
                                                    conditioned_on_parameters) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}
*/

Transform* ProposalKernel::get_transform() const
{
  return this->transform;
}

int ProposalKernel::get_instance_index() const
{
  return instance_index;
}

/*
double ProposalKernel::evaluate_kernel(const Particle &proposed_particle,
                                       const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.parameters,
                                 old_particle.parameters,
                                 this->proposal_parameters);
}

double ProposalKernel::evaluate_kernel(const Particle &proposed_particle,
                                       const Particle &old_particle,
                                       const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.parameters.merge(conditioned_on_parameters),
                                 old_particle.parameters.merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

arma::mat ProposalKernel::gradient_of_log(const std::string &variable,
                                          const Particle &proposed_particle,
                                          const Particle &old_particle)
{
  arma::mat result = this->specific_gradient_of_log(variable,
                                                    proposed_particle,
                                                    old_particle);
  if (this->transform!=NULL)
    result = result * this->transform->jacobian(proposed_particle.parameters);
  
  return result;
}

/*
arma::mat ProposalKernel::gradient_of_log(const std::string &variable,
                                          Particle &proposed_particle,
                                          Particle &old_particle,
                                          const Parameters &conditioned_on_parameters)
{
  arma::mat result = this->specific_gradient_of_log(variable,
                                                    proposed_particle,
                                                    old_particle,
                                                    conditioned_on_parameters);
  if (this->transform!=NULL)
    result = result * this->transform->jacobian(proposed_particle.parameters);

  return result;
}
*/

arma::mat ProposalKernel::subsample_gradient_of_log(const std::string &variable,
                                                    const Particle &proposed_particle,
                                                    const Particle &old_particle)
{
  arma::mat result = this->specific_subsample_gradient_of_log(variable,
                                                              proposed_particle,
                                                              old_particle);
  if (this->transform!=NULL)
    result = result * this->transform->jacobian(proposed_particle.parameters);
  
  return result;
}

/*
arma::mat ProposalKernel::subsample_gradient_of_log(const std::string &variable,
                                                    Particle &proposed_particle,
                                                    Particle &old_particle,
                                                    const Parameters &conditioned_on_parameters)
{
  arma::mat result = this->specific_subsample_gradient_of_log(variable,
                                                              proposed_particle,
                                                              old_particle,
                                                              conditioned_on_parameters);
  if (this->transform!=NULL)
    result = result * this->transform->jacobian(proposed_particle.parameters);
  
  return result;
}
*/
