#include "deterministic_scan_mcmc.h"
#include "deterministic_scan_standard_mcmc_output.h"

DeterministicScanMCMC::DeterministicScanMCMC()
  :MCMC()
{
}

DeterministicScanMCMC::DeterministicScanMCMC(const std::vector<MCMC*> &moves_in)
:MCMC()
{
  this->moves = moves_in;
  this->order.clear();
  this->order.reserve(this->moves.size());
  for (size_t i=0; i<this->moves.size(); ++i)
  {
    this->order.push_back(i);
  }
}

DeterministicScanMCMC::DeterministicScanMCMC(MCMCTermination* termination_in,
                                             const std::vector<MCMC*> &moves_in)
:MCMC(termination_in)
{
  this->moves = moves_in;
  this->order.clear();
  this->order.reserve(this->moves.size());
  for (size_t i=0; i<this->moves.size(); ++i)
  {
    this->order.push_back(i);
  }
}

DeterministicScanMCMC::~DeterministicScanMCMC()
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

//Copy constructor for the DeterministicScanMCMC class.
DeterministicScanMCMC::DeterministicScanMCMC(const DeterministicScanMCMC &another)
  :MCMC(another)
{
  this->make_copy(another);
}

void DeterministicScanMCMC::operator=(const DeterministicScanMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->moves.clear();
  this->order.clear();
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* DeterministicScanMCMC::duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

MCMC* DeterministicScanMCMC::mcmc_duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

DeterministicScanMCMC* DeterministicScanMCMC::deterministic_scan_mcmc_duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

void DeterministicScanMCMC::make_copy(const DeterministicScanMCMC &another)
{
  this->moves.resize(0);
  this->moves.reserve(another.moves.size());
  for (std::vector<MCMC*>::const_iterator i=another.moves.begin();
       i!=another.moves.end();
       ++i)
  {
    if (*i!=NULL)
      this->moves.push_back((*i)->mcmc_duplicate());
    else
      this->moves.push_back(NULL);
  }
  
  this->order = another.order;
}

Particle DeterministicScanMCMC::move(RandomNumberGenerator &rng,
                                     const Particle &particle) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->move(rng,
                                             current_particle);
  }
  return current_particle;
}

/*
Particle DeterministicScanMCMC::move(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     const Parameters &conditioned_on_parameters) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->move(rng,
                                             current_particle,
                                             conditioned_on_parameters);
  }
  return current_particle;
}
*/

Particle DeterministicScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->subsample_move(rng,
                                                       current_particle);
  }
  return current_particle;
}

/*
Particle DeterministicScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     const Parameters &conditioned_on_parameters) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->subsample_move(rng,
                                             current_particle,
                                             conditioned_on_parameters);
  }
  return current_particle;
}
*/

void DeterministicScanMCMC::smc_adapt(SMCOutput* current_state)
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->smc_adapt(current_state);
  }
}

void DeterministicScanMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->ensemble_adapt(current_state);
  }
}

void DeterministicScanMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                                size_t iteration_counter)
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->specific_mcmc_adapt(current_particle,
                              iteration_counter);
  }
}

void DeterministicScanMCMC::set_index(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index(index_in);
  }
}

void DeterministicScanMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}

std::vector<const ProposalKernel*> DeterministicScanMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> all_proposals;
  
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    std::vector<const ProposalKernel*> next_proposals = (*i)->get_proposals();
    if (all_proposals.size()==0)
    {
      all_proposals.insert(all_proposals.end(), next_proposals.begin(), next_proposals.end());
    }
  }
  
  return all_proposals;
}

std::vector<MCMC*> DeterministicScanMCMC::get_duplicate_moves() const
{
  std::vector<MCMC*> duplicate_moves;
  duplicate_moves.reserve(this->moves.size());
  
  for (std::vector<MCMC*>::const_iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    duplicate_moves.push_back((*i)->mcmc_duplicate());
  }
  return duplicate_moves;
}

StandardMCMCOutput* DeterministicScanMCMC::initialise_mcmc_output() const
{
  return new DeterministicScanStandardMCMCOutput(this->deterministic_scan_mcmc_duplicate());
}
