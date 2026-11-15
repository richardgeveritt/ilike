#include "metropolis_hastings_standard_mcmc_output.h"
#include "utils.h"
#include "metropolis_hastings_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput()
:StandardMCMCOutput()
{
  this->mcmc = NULL;
}

MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->termination->set_parameters(this);
}

MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->owns_mcmc = owns_mcmc_in;
  this->termination->set_parameters(this);
}

MetropolisHastingsStandardMCMCOutput::~MetropolisHastingsStandardMCMCOutput()
{
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
}

//Copy constructor for the MetropolisHastingsStandardMCMCOutput class.
MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput(const MetropolisHastingsStandardMCMCOutput &another)
:StandardMCMCOutput(another)
{
  this->make_copy(another);
}

void MetropolisHastingsStandardMCMCOutput::operator=(const MetropolisHastingsStandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* MetropolisHastingsStandardMCMCOutput::duplicate() const
{
  return( new MetropolisHastingsStandardMCMCOutput(*this));
}

void MetropolisHastingsStandardMCMCOutput::make_copy(const MetropolisHastingsStandardMCMCOutput &another)
{
  this->owns_mcmc = another.owns_mcmc;
  if (another.mcmc==NULL)
    this->mcmc = NULL;
  else if (another.owns_mcmc)
    this->mcmc = another.mcmc->metropolis_hastings_mcmc_duplicate();
  else
    this->mcmc = another.mcmc;
}

/*
 void MetropolisHastingsStandardMCMCOutput::specific_mcmc_adapt()
 {
 this->mcmc->mcmc_adapt(this->output.back(),
 this->iteration_counter);
 }
 */

/*
 Particle MetropolisHastingsStandardMCMCOutput::move(RandomNumberGenerator &rng,
 const Particle &particle) const
 {
 return this->mcmc->move(rng, particle);
 }
 
 Particle MetropolisHastingsStandardMCMCOutput::subsample_move(RandomNumberGenerator &rng,
 const Particle &particle) const
 {
 return this->mcmc->subsample_move(rng, particle);
 }
 */

MCMC* MetropolisHastingsStandardMCMCOutput::get_mcmc()
{
  return this->mcmc;
}

const MCMC* MetropolisHastingsStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
}
