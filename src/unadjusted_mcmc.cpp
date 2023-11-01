#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "unadjusted_mcmc.h"
#include "symmetric_proposal_kernel.h"
#include "ensemble_kalman_output.h"
#include "unadjusted_standard_mcmc_output.h"

UnadjustedMCMC::UnadjustedMCMC()
  :MCMC()
{
  this->index = NULL;
  this->proposal = NULL;
}

UnadjustedMCMC::UnadjustedMCMC(size_t number_of_iterations_in,
                               ProposalKernel* proposal_in)
  :MCMC(number_of_iterations_in)
{
  this->proposal = proposal_in;
  this->index = NULL;
}

UnadjustedMCMC::UnadjustedMCMC(size_t number_of_iterations_in,
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

UnadjustedMCMC::UnadjustedMCMC(MCMCTermination* termination_in,
                               ProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal = proposal_in;
}

UnadjustedMCMC::~UnadjustedMCMC()
{
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->index!=NULL)
    delete index;
}

//Copy constructor for the UnadjustedMCMC class.
UnadjustedMCMC::UnadjustedMCMC(const UnadjustedMCMC &another)
  :MCMC(another)
{
  this->make_copy(another);
}

void UnadjustedMCMC::operator=(const UnadjustedMCMC &another)
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

Kernel* UnadjustedMCMC::duplicate() const
{
  return( new UnadjustedMCMC(*this));
}

MCMC* UnadjustedMCMC::mcmc_duplicate() const
{
  return( new UnadjustedMCMC(*this));
}

UnadjustedMCMC* UnadjustedMCMC::unadjusted_mcmc_duplicate() const
{
  return( new UnadjustedMCMC(*this));
}

void UnadjustedMCMC::make_copy(const UnadjustedMCMC &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

Particle UnadjustedMCMC::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  return this->proposal->move(rng,particle);
}

/*
Particle UnadjustedMCMC::move(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.evaluate_likelihoods(this->index,
                                                     conditioned_on_parameters) -
      particle.target_evaluated)
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

Particle UnadjustedMCMC::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  return this->proposal->subsample_move(rng,particle);
}

/*
Particle UnadjustedMCMC::subsample_move(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters) const
{
  Particle proposed_particle = this->proposal->subsample_move(rng,
                                                              particle,
                                                              conditioned_on_parameters);
  
  double log_u = log(runif(rng));
  
  if (log_u < proposed_particle.subsample_evaluate_likelihoods(this->index,
                                                               conditioned_on_parameters) -
      particle.subsample_target_evaluated)
  {
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}
*/

void UnadjustedMCMC::smc_adapt(SMCOutput* current_state)
{
  proposal->smc_adapt(current_state);
}

void UnadjustedMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  this->proposal->ensemble_adapt(current_state);
}

void UnadjustedMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                                 size_t iteration_counter)
{
  this->proposal->mcmc_adapt(current_particle,
                             iteration_counter);
}

void UnadjustedMCMC::set_index(Index* index_in)
{
  if (this->index!=NULL)
    delete this->index;
  
  this->index = index_in;
  this->proposal->set_index(index_in);
}

void UnadjustedMCMC::set_index_if_null(Index* index_in)
{
  if (this->index==NULL)
    this->index = index_in;
  this->proposal->set_index_if_null(index_in);
}

void UnadjustedMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

std::vector<const ProposalKernel*> UnadjustedMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals = this->proposal->get_proposals();
  //proposals.push_back(this->proposal);
  return proposals;
}

StandardMCMCOutput* UnadjustedMCMC::initialise_mcmc_output() const
{
  return new UnadjustedStandardMCMCOutput(this->unadjusted_mcmc_duplicate());
}
