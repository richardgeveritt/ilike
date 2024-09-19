#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "metropolis_hastings_mcmc.h"
#include "gaussian_random_walk_proposal_kernel.h"
#include "metropolis_hastings_standard_mcmc_output.h"

namespace ilike
{
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
                                               const std::string &variable_name_in,
                                               const arma::mat &proposal_covariance_in)
:MCMC(number_of_iterations_in)
{
  this->index = NULL;
  std::vector<std::string> variable_names_in;
  variable_names_in.push_back(variable_name_in);
  std::vector<arma::mat> proposal_covariances_in;
  proposal_covariances_in.push_back(proposal_covariance_in);
  this->proposal = new GaussianRandomWalkProposalKernel(variable_names_in,
                                                        proposal_covariances_in);
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               const std::vector<std::string> &variable_names_in,
                                               const std::vector<arma::mat> &proposal_covariances_in)
:MCMC(number_of_iterations_in)
{
  this->index = NULL;
  this->proposal = new GaussianRandomWalkProposalKernel(variable_names_in,
                                                        proposal_covariances_in);
  // default to Gaussian random walk
  //this->proposal = ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in);
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(MCMCTermination* termination_in,
                                               ProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal = proposal_in;
  this->index = NULL;
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

MetropolisHastingsMCMC* MetropolisHastingsMCMC::metropolis_hastings_mcmc_duplicate() const
{
  return( new MetropolisHastingsMCMC(*this));
}

void MetropolisHastingsMCMC::make_copy(const MetropolisHastingsMCMC &another)
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

Particle MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
                                      const Particle &particle) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle);
  
  double log_u = log(runif(rng));
  
  double numerator = proposed_particle.evaluate_likelihoods(this->index);
  double denominator = particle.target_evaluated;
  
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
  {
    return particle;
  }
  
  if (log_u < numerator - denominator +
      this->proposal->evaluate_kernel(particle, proposed_particle) -
      this->proposal->evaluate_kernel(proposed_particle, particle))
  {
    proposed_particle.set_acceptance(this->proposal,true);
    return proposed_particle;
  }
  else
  {
    //particle.set_acceptance(this->proposal,false);
    return particle;
  }
  
}

/*
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
 */

Particle MetropolisHastingsMCMC::subsample_move(RandomNumberGenerator &rng,
                                                const Particle &particle) const
{
  Particle proposed_particle = this->proposal->subsample_move(rng,
                                                              particle);
  
  double log_u = log(runif(rng));
  
  double numerator = proposed_particle.subsample_evaluate_likelihoods(this->index);
  double denominator = particle.subsample_target_evaluated;
  
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
  {
    return particle;
  }
  
  if (log_u < numerator - denominator +
      this->proposal->subsample_evaluate_kernel(particle, proposed_particle) -
      this->proposal->subsample_evaluate_kernel(proposed_particle, particle))
  {
    proposed_particle.set_acceptance(this->proposal,true);
    return proposed_particle;
  }
  else
  {
    //particle.set_acceptance(this->proposal,false);
    return particle;
  }
  
}

/*
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
 */

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

void MetropolisHastingsMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                                 size_t iteration_counter)
{
  this->proposal->mcmc_adapt(current_particle,
                             iteration_counter);
}

void MetropolisHastingsMCMC::set_index(Index* index_in)
{
  if (this->index!=NULL)
    delete this->index;
  
  this->index = index_in;
  this->proposal->set_index(index_in);
}

void MetropolisHastingsMCMC::set_index_if_null(Index* index_in)
{
  if (this->index==NULL)
    this->index = index_in;
  if (this->proposal!=NULL)
  {
    this->proposal->set_index_if_null(index_in);
  }
}

void MetropolisHastingsMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  if (this->proposal!=NULL)
  {
    this->proposal->set_proposal_parameters(proposal_parameters_in);
  }
}

std::vector<const ProposalKernel*> MetropolisHastingsMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals;
  if (this->proposal!=NULL)
  {
    proposals = this->proposal->get_proposals();
  }
  //proposals.push_back(this->proposal);
  return proposals;
}

StandardMCMCOutput* MetropolisHastingsMCMC::initialise_mcmc_output() const
{
  return new MetropolisHastingsStandardMCMCOutput(this->metropolis_hastings_mcmc_duplicate());
}
}
