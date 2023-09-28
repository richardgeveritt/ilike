#include "stochastic_scan_mcmc.h"
#include "distributions.h"
#include "stochastic_scan_standard_mcmc_output.h"

StochasticScanMCMC::StochasticScanMCMC()
  :MCMC()
{
}

StochasticScanMCMC::StochasticScanMCMC(const std::vector<MCMC*> &moves_in,
                                       const arma::colvec &unnormalised_probabilities_in)
:MCMC()
{
  this->moves = moves_in;
  this->probabilities = unnormalised_probabilities_in/sum(unnormalised_probabilities_in);
}

StochasticScanMCMC::StochasticScanMCMC(MCMCTermination* termination_in,
                                       const std::vector<MCMC*> &moves_in,
                                       const arma::colvec &unnormalised_probabilities_in)
:MCMC(termination_in)
{
  this->moves = moves_in;
  this->probabilities = unnormalised_probabilities_in/sum(unnormalised_probabilities_in);
}

StochasticScanMCMC::~StochasticScanMCMC()
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

//Copy constructor for the StochasticScanMCMC class.
StochasticScanMCMC::StochasticScanMCMC(const StochasticScanMCMC &another)
  :MCMC(another)
{
  this->make_copy(another);
}

void StochasticScanMCMC::operator=(const StochasticScanMCMC &another)
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
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* StochasticScanMCMC::duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

MCMC* StochasticScanMCMC::mcmc_duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

StochasticScanMCMC* StochasticScanMCMC::stochastic_scan_mcmc_duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

void StochasticScanMCMC::make_copy(const StochasticScanMCMC &another)
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
  
  this->probabilities = another.probabilities;
}

Particle StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                  const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle);
}

/*
Particle StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                  Particle &particle,
                                  const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle,
                                                           conditioned_on_parameters);
}
*/

Particle StochasticScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                            const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     particle);
}

/*
Particle StochasticScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                            Particle &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     particle,
                                                                     conditioned_on_parameters);
}
*/

void StochasticScanMCMC::smc_adapt(SMCOutput* current_state)
{
  // Adapt probabilities?
  
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->smc_adapt(current_state);
  }
}

/*
EnsembleMember StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                        const Index* index,
                                        EnsembleMember &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           index,
                                                           particle);
}

EnsembleMember StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                        const Index* index,
                                        EnsembleMember &particle,
                                        const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           index,
                                                           particle,
                                                           conditioned_on_parameters);
}

EnsembleMember StochasticScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                                  const Index* index,
                                                  EnsembleMember &particle,
                                                  const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     index,
                                                                     particle,
                                                                     conditioned_on_parameters);
}
*/
 
void StochasticScanMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  // Adapt probabilities?
  
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->ensemble_adapt(current_state);
  }
}

void StochasticScanMCMC::specific_mcmc_adapt(const Particle &current_particle,
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

void StochasticScanMCMC::set_index(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index(index_in);
  }
}

void StochasticScanMCMC::set_index_if_null(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index_if_null(index_in);
  }
}

void StochasticScanMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}

std::vector<const ProposalKernel*> StochasticScanMCMC::get_proposals() const
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

std::vector<MCMC*> StochasticScanMCMC::get_duplicate_moves() const
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

StandardMCMCOutput* StochasticScanMCMC::initialise_mcmc_output() const
{
  return new StochasticScanStandardMCMCOutput(this->stochastic_scan_mcmc_duplicate());
}
