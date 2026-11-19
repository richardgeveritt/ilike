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

namespace ilike
{
int ProposalKernel::instance_counter = 0;

ProposalKernel::ProposalKernel()
:Kernel()
{
  this->instance_index = this->instance_counter++;
  this->mcmc_adaptor.reset();
  this->smc_adaptor.reset();
  this->transform.reset();
}

ProposalKernel::~ProposalKernel()
{
  this->mcmc_adaptor.reset();
  
  this->smc_adaptor.reset();
  
  this->transform.reset();
}

ProposalKernel::ProposalKernel(const Parameters &proposal_parameters_in)
:Kernel()
{
  this->instance_index = this->instance_counter++;
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
  
  this->mcmc_adaptor.reset();
  
  this->smc_adaptor.reset();
  
  this->transform.reset();
  
  Kernel::operator=(another);
  this->make_copy(another);
}

void ProposalKernel::make_copy(const ProposalKernel &another)
{
  if (another.mcmc_adaptor!=NULL)
  {
    // the copy's adaptor adapts the copy, not the kernel copied from
    this->mcmc_adaptor.reset(another.mcmc_adaptor->duplicate());
    this->mcmc_adaptor->set_proposal(this);
  }
  else
    this->mcmc_adaptor.reset();
  
  if (another.smc_adaptor!=NULL)
  {
    // the copy's adaptor adapts the copy, not the kernel copied from
    this->smc_adaptor.reset(another.smc_adaptor->duplicate());
    this->smc_adaptor->set_proposal(this);
  }
  else
    this->smc_adaptor.reset();
  
  if (another.transform!=NULL)
    this->transform.reset(another.transform->duplicate());
  else
    this->transform.reset();

  this->instance_index = another.instance_index;
  this->instance_counter = another.instance_counter;
}

Particle ProposalKernel::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  return this->move(rng, particle, NULL);
}

Particle ProposalKernel::move(RandomNumberGenerator &rng,
                              const Particle &particle,
                              const Index* new_outputs_index) const
{
  return this->proposed_from(particle, this->simulate(rng,particle), new_outputs_index);
}

// the proposed particle for parameters simulated by the kernel (on its own,
// possibly transformed, space): the steps of move() after the simulation
Particle ProposalKernel::proposed_from(const Particle &particle,
                                       const Parameters &simulated,
                                       const Index* new_outputs_index) const
{
  Particle proposed_particle;
  proposed_particle.parameters = particle.parameters;
  
  if (this->transform==NULL)
  {
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(simulated);
  }
  else
  {
    proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform->inverse_transform(simulated));
  }
  
  // Outputs are created here, with memory managed by Particle hereafter.
  if (particle.factor_variables!=NULL)
  {
    if (new_outputs_index!=NULL)
    {
      proposed_particle.simulate_factor_variables(particle.factor_variables.get(),
                                                  new_outputs_index);
    }
    else
    {
      const Factors* old_factors = particle.factor_variables->get_factors();
      if (old_factors!=NULL)
        proposed_particle.simulate_factor_variables(old_factors);
    }
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

void ProposalKernel::move_into(RandomNumberGenerator &rng,
                               const Particle &particle,
                               const Index* new_outputs_index,
                               Particle &scratch) const
{
  // the steps of move(), in place: the parameters are the current ones with
  // the proposed variables overwritten, written into scratch's own matrices
  // (ProposalTarget; simulate reads only the current particle, so the result
  // and the random numbers are those of move())
  ProposalTarget target(scratch.parameters, particle.parameters);
  if (this->transform==NULL)
  {
    this->simulate_into(rng, particle, target);
  }
  else
  {
    Parameters proposed_parameters_in_transformed_space = this->simulate(rng,particle);
    target.write_all(this->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  target.finish();
  this->finish_into(particle, new_outputs_index, scratch);
}

// move_into for parameters simulated by the kernel (on its own space)
void ProposalKernel::into_from(const Particle &particle,
                               const Parameters &simulated,
                               const Index* new_outputs_index,
                               Particle &scratch) const
{
  ProposalTarget target(scratch.parameters, particle.parameters);
  if (this->transform==NULL)
    target.write_all(simulated);
  else
    target.write_all(this->transform->inverse_transform(simulated));
  target.finish();
  this->finish_into(particle, new_outputs_index, scratch);
}

// the steps of move_into after the parameters are written
void ProposalKernel::finish_into(const Particle &particle,
                                 const Index* new_outputs_index,
                                 Particle &scratch) const
{
  bool recycled = false;
  if ((particle.factor_variables!=NULL) && (scratch.factor_variables!=NULL))
  {
    if (new_outputs_index!=NULL)
      recycled = scratch.factor_variables->recycle_with_new_outputs(particle.factor_variables.get(),
                                                                    scratch.parameters,
                                                                    new_outputs_index);
    else if (particle.factor_variables->get_factors()!=NULL)
      recycled = scratch.factor_variables->recycle_all_outputs(particle.factor_variables->get_factors(),
                                                               scratch.parameters);
    if (recycled)
      scratch.factor_variables->set_particle(&scratch);
  }
  if (!recycled)
  {
    if (scratch.factor_variables!=NULL)
    {
      scratch.factor_variables.reset();
      scratch.factor_variables.reset();
    }
    if (particle.factor_variables!=NULL)
    {
      if (new_outputs_index!=NULL)
      {
        scratch.simulate_factor_variables(particle.factor_variables.get(),
                                          new_outputs_index);
      }
      else
      {
        const Factors* old_factors = particle.factor_variables->get_factors();
        if (old_factors!=NULL)
          scratch.simulate_factor_variables(old_factors);
      }
    }
  }

  if (scratch.ensemble_factor_variables!=NULL)
  {
    scratch.ensemble_factor_variables.reset();
    scratch.ensemble_factor_variables.reset();
  }
  if (particle.ensemble_factor_variables!=NULL)
  {
    const EnsembleFactors* old_ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (old_ensemble_factors!=NULL)
      scratch.simulate_ensemble_factor_variables(old_ensemble_factors);
  }

  // as a new Particle: empty proposal store, targets not yet evaluated
  scratch.current_proposal_store.clear();
  scratch.previous_target_evaluated = 0.0;
  scratch.target_evaluated = 0.0;
  scratch.subsample_previous_target_evaluated = 0.0;
  scratch.subsample_target_evaluated = 0.0;
  scratch.previous_ensemble_target_evaluated = 0.0;
  scratch.subsample_previous_ensemble_target_evaluated = 0.0;
  scratch.ensemble_target_evaluated = 0.0;
  scratch.subsample_ensemble_target_evaluated = 0.0;
  scratch.simulate_proposal_variables(particle.proposals_to_transform_for_pointer, particle.proposals_to_find_gradient_for_pointer);
  scratch.accepted_outputs = particle.accepted_outputs;
  scratch.previous_self = &particle;
}

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

  return proposed_particle; // what we need is for previous target eval to be set to eval for this target (when evaluated there)!!!!! we don't want it to know its grad, until evaluated there
}

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

void ProposalKernel::set_smc_adaptor(SMCAdaptor* adaptor_in)
{
  this->smc_adaptor.reset(adaptor_in);
  if (this->smc_adaptor!=NULL)
    this->smc_adaptor->set_proposal(this);
}

void ProposalKernel::set_mcmc_adaptor(MCMCAdaptor* adaptor_in)
{
  this->mcmc_adaptor.reset(adaptor_in);
  if (this->mcmc_adaptor!=NULL)
    this->mcmc_adaptor->set_proposal(this);
}

void ProposalKernel::mcmc_adapt(const AdaptationStep &step)
{
  if (this->mcmc_adaptor!=NULL)
    this->mcmc_adaptor->adapt_unless_stopped(step);
}

void ProposalKernel::prepare_mcmc_adaptation()
{
  if (this->mcmc_adaptor!=NULL)
    this->mcmc_adaptor->prepare_for_move();
}

void ProposalKernel::append_adaptation_state(std::vector<std::string> *names,
                                             std::vector<double> *values,
                                             bool final,
                                             const std::string &prefix) const
{
  if (this->mcmc_adaptor!=NULL)
  {
    this->mcmc_adaptor->append_state(names, values, final, prefix);
    // with until =: the adaptation step at which it stopped (0: it did not)
    if (final && this->mcmc_adaptor->has_until())
    {
      if (names!=NULL)
        names->push_back(prefix + this->mcmc_adaptor->state_label() + ".adaptation_stopped_at");
      if (values!=NULL)
        values->push_back(double(this->mcmc_adaptor->stopped_at()));
    }
  }
}

Parameters* ProposalKernel::get_proposal_parameters() const
{
  return NULL;
}

bool ProposalKernel::proposal_uses_gaussian_mean() const
{
  return false;
}

GaussianProposalInfo* ProposalKernel::get_gaussian_proposal_info(const std::string &variable)
{
  return NULL;
}

arma::colvec* ProposalKernel::get_component_scales(const std::string &variable)
{
  return NULL;
}

PrincipalComponentsInfo* ProposalKernel::get_principal_components(const std::string &variable)
{
  return NULL;
}

void ProposalKernel::mcmc_adapt(const Particle &current_particle,
                                size_t iteration_counter)
{
  if (this->mcmc_adaptor!=NULL)
    this->mcmc_adaptor->mcmc_adapt(current_particle,
                                   iteration_counter);
}

double ProposalKernel::evaluate_kernel(const Particle &proposed_particle,
                                       const Particle &old_particle) const
{
  if (this->transform==NULL)
  {
    return this->specific_evaluate_kernel(proposed_particle, old_particle);
  }
  else
  {
    
    // need to include absolute value of Jacobian determinant
    return this->specific_evaluate_kernel(proposed_particle, old_particle) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}

void ProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                   const Particle &particle,
                                   ProposalTarget &target) const
{
  target.write_all(this->simulate(rng, particle));
}

double ProposalKernel::evaluate_kernel(const Parameters &proposed_parameters,
                                       const Parameters &old_parameters) const
{
  if (this->transform==NULL)
    return this->specific_evaluate_kernel_from_parameters(proposed_parameters, old_parameters);
  else
    return this->evaluate_kernel(Particle(proposed_parameters), Particle(old_parameters));
}

double ProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                const Parameters &old_parameters) const
{
  return this->specific_evaluate_kernel(Particle(proposed_parameters), Particle(old_parameters));
}

double ProposalKernel::subsample_evaluate_kernel(const Particle &proposed_particle,
                                                 const Particle &old_particle) const
{
  if (this->transform==NULL)
  {
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle);
  }
  else
  {
    
    // need to include absolute value of Jacobian determinant
    return this->specific_subsample_evaluate_kernel(proposed_particle,
                                                    old_particle) + this->transform->log_abs_jacobian_determinant(proposed_particle.parameters);
  }
}

Transform* ProposalKernel::get_transform() const
{
  return this->transform.get();
}

int ProposalKernel::get_instance_index() const
{
  return instance_index;
}

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

}
