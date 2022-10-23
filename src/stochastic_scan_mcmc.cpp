#include "stochastic_scan_mcmc.h"
#include "distributions.h"

StochasticScanMCMC::StochasticScanMCMC()
  :MCMC()
{
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
                                  Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle);
}

Particle StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                  Particle &particle,
                                  const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle,
                                                           conditioned_on_parameters);
}

Particle StochasticScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                            Particle &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     particle,
                                                                     conditioned_on_parameters);
}

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

void StochasticScanMCMC::specific_mcmc_adapt(Particle &current_particle,
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
