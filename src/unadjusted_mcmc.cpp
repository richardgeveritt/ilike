#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "unadjusted_mcmc.h"
#include "symmetric_proposal_kernel.h"
#include "ensemble_kalman_output.h"
#include "unadjusted_standard_mcmc_output.h"

namespace ilike
{
UnadjustedMCMC::UnadjustedMCMC()
:MCMC()
{
  this->index.reset();
  this->proposal.reset();
}

UnadjustedMCMC::UnadjustedMCMC(size_t number_of_iterations_in,
                               ProposalKernel* proposal_in)
:MCMC(number_of_iterations_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

UnadjustedMCMC::UnadjustedMCMC(size_t number_of_iterations_in,
                               const std::vector<Parameters> &initial_points_in,
                               const Parameters &proposal_variances_in)
:MCMC(number_of_iterations_in)
{
  this->index.reset();
  this->proposal.reset();
  // default to Gaussian random walk
  //this->proposal.reset(ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in));
}

UnadjustedMCMC::UnadjustedMCMC(MCMCTermination* termination_in,
                               ProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

UnadjustedMCMC::~UnadjustedMCMC()
{
  this->proposal.reset();
  
  this->index.reset();
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
  
  this->proposal.reset();
  
  this->index.reset();
  
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
    this->proposal.reset(another.proposal->proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

Particle UnadjustedMCMC::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  // No factor is evaluated: avoiding that is often the reason for using an
  // unadjusted move, and the user may not supply the functions to do it. The
  // moved particle gets new (unevaluated) outputs for this move's factors and
  // cut factors, and carries the others over. make_mcmc refuses to combine
  // this move with a move whose acceptance ratio would read the unevaluated
  // values (check_unadjusted_compatibility).
  return this->proposal->move(rng,
                              particle,
                              this->get_new_outputs_index(this->index.get()));
}

Particle UnadjustedMCMC::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  return this->proposal->subsample_move(rng,particle);
}

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

// Ownership rule for indices: set_index and set_index_if_null never take
// ownership of index_in. This object keeps its own duplicate (deleted in the
// destructor), and the caller remains responsible for index_in. Previously the
// pointer itself was stored, so when one index was handed to several moves (e.g.
// the moves of a sweep) each destructor deleted it -> double free.
void UnadjustedMCMC::set_index(Index* index_in)
{
  this->index.reset();
  this->index.reset((index_in!=NULL) ? index_in->duplicate() : NULL);
  if (this->proposal!=NULL)
    this->proposal->set_index(index_in);
  this->update_index_with_cut(this->index.get());
}

void UnadjustedMCMC::set_index_if_null(Index* index_in)
{
  if ( (this->index==NULL) && (index_in!=NULL) )
    this->index.reset(index_in->duplicate());
  if (this->proposal!=NULL)
    this->proposal->set_index_if_null(index_in);
  this->update_index_with_cut(this->index.get());
}

bool UnadjustedMCMC::reads_cached_factors() const
{
  return false;
}

void UnadjustedMCMC::set_cut_index(const Index* cut_index_in)
{
  this->store_cut_index(cut_index_in);
  this->update_index_with_cut(this->index.get());
}

void UnadjustedMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

std::vector<const ProposalKernel*> UnadjustedMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals = this->proposal->get_proposals();
  return proposals;
}

bool UnadjustedMCMC::adapts() const
{
  return this->proposal->has_mcmc_adaptor();
}

StandardMCMCOutput* UnadjustedMCMC::initialise_mcmc_output() const
{
  // A move that adapts gives each chain its own copy of the kernel; otherwise
  // all the particles' outputs share this one (used read-only by the moves).
  if (this->adapts())
    return new UnadjustedStandardMCMCOutput(this->unadjusted_mcmc_duplicate());
  return new UnadjustedStandardMCMCOutput(const_cast<UnadjustedMCMC*>(this), false);
}
}
