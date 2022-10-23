#include <iostream>
#include "parameter_particle_simulator.h"
#include "likelihood_estimator.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"
#include "distributions.h"
#include "independent_proposal_kernel.h"
#include "custom_independent_proposal_kernel.h"
#include "factors.h"
#include "vector_factor_variables.h"

//Default constructor.
ParameterParticleSimulator::ParameterParticleSimulator()
  :ParticleSimulator()
{
  this->proposal = NULL;
}

/*
ParameterParticleSimulator::ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                                                       const std::string &resample_variable_name_in)
  :ParticleSimulator(resample_variable_name_in)
{
  this->simulate_parameters = simulate_parameters_in;
  this->likelihood_estimators = likelihood_estimators_in;
}
*/

ParameterParticleSimulator::ParameterParticleSimulator(IndependentProposalKernel* proposal_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:ParticleSimulator()
{
  this->proposal = proposal_in;
  this->likelihood_estimators = likelihood_estimators_in;
}

//Copy constructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::ParameterParticleSimulator(const ParameterParticleSimulator &another)
  :ParticleSimulator(another)
{
  this->make_copy(another);
}

//Destructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::~ParameterParticleSimulator()
{
  if (this->proposal!=NULL)
    delete this->proposal;
}

void ParameterParticleSimulator::operator=(const ParameterParticleSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  this->likelihood_estimators.clear();

  ParticleSimulator::operator=(another);
  this->make_copy(another);
}

ParticleSimulator* ParameterParticleSimulator::duplicate()const
{
  return( new ParameterParticleSimulator(*this));
}

void ParameterParticleSimulator::make_copy(const ParameterParticleSimulator &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->independent_proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  this->likelihood_estimators = another.likelihood_estimators;
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              Factors* factors) const
{
  Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //std::cout<<simulated_parameters<<std::endl;
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, simulated_factor_variables);
  
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              Factors* factors,
                                              const Parameters &conditioned_on_parameters) const
{
  Parameters simulated_parameters = this->proposal->independent_simulate(rng,
                                                                         conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters,
                                                                                   conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, simulated_factor_variables);
  
}

Particle ParameterParticleSimulator::subsample_simulate(RandomNumberGenerator &rng,
                                                        Factors* factors,
                                                        const Parameters &conditioned_on_parameters) const
{
  Parameters simulated_parameters = this->proposal->independent_simulate(rng,
                                                                         conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  FactorVariables* simulated_factor_variables = factors->subsample_simulate_factor_variables(simulated_parameters,
                                                                                             conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, simulated_factor_variables);
  
}

double ParameterParticleSimulator::evaluate(Particle &input) const
{
  // could generalise to evaluating the proposals for the likelihoodestimators
  return this->proposal->evaluate_independent_kernel(input.parameters);
}

double ParameterParticleSimulator::evaluate(Particle &input,
                                            const Parameters &conditioned_on_parameters) const
{
    // could generalise to evaluating the proposals for the likelihoodestimators
  return this->proposal->evaluate_independent_kernel(input.parameters.merge(conditioned_on_parameters));
}

double ParameterParticleSimulator::subsample_evaluate(Particle &input,
                                                      const Parameters &conditioned_on_parameters) const
{
  // could generalise to evaluating the proposals for the likelihoodestimators
  return this->proposal->subsample_evaluate_independent_kernel(input.parameters.merge(conditioned_on_parameters));
}
