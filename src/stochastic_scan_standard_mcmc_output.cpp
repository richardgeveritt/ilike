#include "stochastic_scan_standard_mcmc_output.h"
#include "utils.h"
#include "stochastic_scan_mcmc.h"
#include "mcmc_termination.h"

StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput()
  :StandardMCMCOutput()
{
  this->mcmc = NULL;
}

StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput(StochasticScanMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->termination->set_parameters(this);
}

StochasticScanStandardMCMCOutput::~StochasticScanStandardMCMCOutput()
{
  if (this->mcmc!=NULL)
    delete this->mcmc;
}

//Copy constructor for the StochasticScanStandardMCMCOutput class.
StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput(const StochasticScanStandardMCMCOutput &another)
  :StandardMCMCOutput(another)
{
  this->make_copy(another);
}

void StochasticScanStandardMCMCOutput::operator=(const StochasticScanStandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* StochasticScanStandardMCMCOutput::duplicate() const
{
  return( new StochasticScanStandardMCMCOutput(*this));
}

void StochasticScanStandardMCMCOutput::make_copy(const StochasticScanStandardMCMCOutput &another)
{
  if (another.mcmc!=NULL)
    this->mcmc = another.mcmc->stochastic_scan_mcmc_duplicate();
  else
    this->mcmc = NULL;
}

/*
void StochasticScanStandardMCMCOutput::specific_mcmc_adapt()
{
  for (std::vector<MCMC*>::iterator i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->mcmc_adapt(this->output.back(),
                     this->iteration_counter);
  }
}
*/

/*
Particle StochasticScanStandardMCMCOutput::move(RandomNumberGenerator &rng,
                                                const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle);
}

Particle StochasticScanStandardMCMCOutput::subsample_move(RandomNumberGenerator &rng,
                                                          const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     particle);
}
*/

MCMC* StochasticScanStandardMCMCOutput::get_mcmc()
{
  return this->mcmc;
}

const MCMC* StochasticScanStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
