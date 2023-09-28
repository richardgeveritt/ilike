#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "metropolis_mcmc.h"
#include "symmetric_proposal_kernel.h"
#include "ensemble_kalman_output.h"
#include "metropolis_standard_mcmc_output.h"

MetropolisMCMC::MetropolisMCMC()
  :MCMC()
{
  this->index = NULL;
  this->proposal = NULL;
}

MetropolisMCMC::MetropolisMCMC(size_t number_of_iterations_in,
                               SymmetricProposalKernel* proposal_in)
  :MCMC(number_of_iterations_in)
{
  this->proposal = proposal_in;
  this->index = NULL;
}

MetropolisMCMC::MetropolisMCMC(size_t number_of_iterations_in,
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

MetropolisMCMC::MetropolisMCMC(MCMCTermination* termination_in,
                               SymmetricProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal = proposal_in;
}

MetropolisMCMC::~MetropolisMCMC()
{
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->index!=NULL)
    delete index;
}

//Copy constructor for the MetropolisMCMC class.
MetropolisMCMC::MetropolisMCMC(const MetropolisMCMC &another)
  :MCMC(another)
{
  this->make_copy(another);
}

void MetropolisMCMC::operator=(const MetropolisMCMC &another)
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

Kernel* MetropolisMCMC::duplicate() const
{
  return( new MetropolisMCMC(*this));
}

MCMC* MetropolisMCMC::mcmc_duplicate() const
{
  return( new MetropolisMCMC(*this));
}

MetropolisMCMC* MetropolisMCMC::metropolis_mcmc_duplicate() const
{
  return( new MetropolisMCMC(*this));
}

void MetropolisMCMC::make_copy(const MetropolisMCMC &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->symmetric_proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

Particle MetropolisMCMC::move(RandomNumberGenerator &rng,
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
  
  if (log_u < numerator - denominator)
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
Particle MetropolisMCMC::move(RandomNumberGenerator &rng,
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

Particle MetropolisMCMC::subsample_move(RandomNumberGenerator &rng,
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
  
  if (log_u < numerator - denominator)
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
Particle MetropolisMCMC::subsample_move(RandomNumberGenerator &rng,
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

void MetropolisMCMC::smc_adapt(SMCOutput* current_state)
{
  proposal->smc_adapt(current_state);
}

void MetropolisMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  this->proposal->ensemble_adapt(current_state);
}

void MetropolisMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                                 size_t iteration_counter)
{
  this->proposal->mcmc_adapt(current_particle,
                             iteration_counter);
}

void MetropolisMCMC::set_index(Index* index_in)
{
  if (this->index!=NULL)
    delete this->index;
  
  this->index = index_in;
  this->proposal->set_index(index_in);
}

void MetropolisMCMC::set_index_if_null(Index* index_in)
{
  if (this->index==NULL)
    this->index = index_in;
  this->proposal->set_index_if_null(index_in);
}

void MetropolisMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

std::vector<const ProposalKernel*> MetropolisMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals = this->proposal->get_proposals();
  //proposals.push_back(this->proposal);
  return proposals;
}

StandardMCMCOutput* MetropolisMCMC::initialise_mcmc_output() const
{
  return new MetropolisStandardMCMCOutput(this->metropolis_mcmc_duplicate());
}
