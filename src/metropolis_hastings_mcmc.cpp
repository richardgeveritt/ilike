#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "metropolis_hastings_mcmc.h"

MetropolisHastingsMCMC::MetropolisHastingsMCMC()
  :MCMC()
{
  this->index = NULL;
  this->proposal = NULL;
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               ProposalKernel* proposal_in)
  :MCMC(number_of_iterations_in)
{
  this->proposal = proposal_in;
  this->index = NULL;
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               const std::vector<Parameters> &initial_points_in,
                                               const Parameters &proposal_variances_in)
  :MCMC(number_of_iterations_in)
{
  this->index = NULL;
  this->proposal = NULL;
  // default to Gaussian random walk
  //this->proposal = ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in);
}

MetropolisHastingsMCMC::~MetropolisHastingsMCMC()
{
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->index!=NULL)
    delete index;
}

//Copy constructor for the MetropolisHastingsMCMC class.
MetropolisHastingsMCMC::MetropolisHastingsMCMC(const MetropolisHastingsMCMC &another)
  :MCMC(another)
{
  this->make_copy(another);
}

void MetropolisHastingsMCMC::operator=(const MetropolisHastingsMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->index!=NULL)
    delete index;
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* MetropolisHastingsMCMC::duplicate() const
{
  return( new MetropolisHastingsMCMC(*this));
}

MCMC* MetropolisHastingsMCMC::mcmc_duplicate() const
{
  return( new MetropolisHastingsMCMC(*this));
}

void MetropolisHastingsMCMC::make_copy(const MetropolisHastingsMCMC &another)
{
  this->proposal = another.proposal;
  if (another.proposal!=NULL)
    this->proposal = another.proposal->proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

Particle MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
                                      Particle &particle) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.evaluate_likelihoods(this->index) -
      particle.target_evaluated +
      this->proposal->evaluate_kernel(particle, proposed_particle) -
      this->proposal->evaluate_kernel(proposed_particle, particle))
  {
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

Particle MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
                                      Particle &particle,
                                      const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.evaluate_likelihoods(this->index,
                                                     conditioned_on_parameters) -
      particle.target_evaluated +
      this->proposal->evaluate_kernel(particle, proposed_particle, conditioned_on_parameters) -
      this->proposal->evaluate_kernel(proposed_particle, particle, conditioned_on_parameters))
  {
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

Particle MetropolisHastingsMCMC::subsample_move(RandomNumberGenerator &rng,
                                                Particle &particle,
                                                const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle = this->proposal->subsample_move(rng,
                                                              particle,
                                                              conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.subsample_evaluate_likelihoods(this->index,
                                                               conditioned_on_parameters) -
      particle.subsample_target_evaluated +
      this->proposal->subsample_evaluate_kernel(particle, proposed_particle, conditioned_on_parameters) -
      this->proposal->subsample_evaluate_kernel(proposed_particle, particle, conditioned_on_parameters))
  {
    proposed_particle.set_acceptance(this->proposal,true);
    return proposed_particle;
  }
  else
  {
    proposed_particle.set_acceptance(this->proposal,false);
    return particle;
  }
  
}

/*
EnsembleMember MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
                                            const Index* index,
                                            EnsembleMember &particle) const
{
  EnsembleMember proposed_particle = this->proposal->move(rng,
                                                          particle);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.evaluate_all_likelihoods(index) -
      particle.target_evaluated +
      this->proposal->evaluate_kernel(&particle, &proposed_particle) -
      this->proposal->evaluate_kernel(&proposed_particle, &particle))
  {
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

EnsembleMember MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
                                            const Index* index,
                                            EnsembleMember &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  EnsembleMember proposed_particle = this->proposal->move(rng,
                                                          particle,
                                                          conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.evaluate_all_likelihoods(index,
                                                         conditioned_on_parameters) -
      particle.target_evaluated +
      this->proposal->evaluate_kernel(&particle, &proposed_particle, conditioned_on_parameters) -
      this->proposal->evaluate_kernel(&proposed_particle, &particle, conditioned_on_parameters))
  {
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

EnsembleMember MetropolisHastingsMCMC::subsample_move(RandomNumberGenerator &rng,
                                                      const Index* index,
                                                      EnsembleMember &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  EnsembleMember proposed_particle = this->proposal->subsample_move(rng,
                                                              particle,
                                                              conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.subsample_evaluate_all_likelihoods(index,
                                                                   conditioned_on_parameters) -
      particle.subsample_target_evaluated +
      this->proposal->subsample_evaluate_kernel(&particle, &proposed_particle, conditioned_on_parameters) -
      this->proposal->subsample_evaluate_kernel(&proposed_particle, &particle, conditioned_on_parameters))
  {
    proposed_particle.set_acceptance(this->proposal,true);
    return proposed_particle;
  }
  else
  {
    proposed_particle.set_acceptance(this->proposal,false);
    return particle;
  }
  
}
*/

void MetropolisHastingsMCMC::smc_adapt(SMCOutput* current_state)
{
  this->proposal->smc_adapt(current_state);
}

void MetropolisHastingsMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  this->proposal->ensemble_adapt(current_state);
}

void MetropolisHastingsMCMC::specific_mcmc_adapt(Particle &current_particle,
                                                 size_t iteration_counter)
{
  this->proposal->mcmc_adapt(current_particle,
                             iteration_counter);
}
