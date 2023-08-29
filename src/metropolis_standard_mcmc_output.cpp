#include "metropolis_standard_mcmc_output.h"
#include "utils.h"
#include "metropolis_mcmc.h"
#include "mcmc_termination.h"

MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput()
  :StandardMCMCOutput()
{
  this->mcmc = NULL;
}

MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->termination->set_parameters(this);
}

MetropolisStandardMCMCOutput::~MetropolisStandardMCMCOutput()
{
  if (this->mcmc!=NULL)
    delete this->mcmc;
}

//Copy constructor for the MetropolisStandardMCMCOutput class.
MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput(const MetropolisStandardMCMCOutput &another)
  :StandardMCMCOutput(another)
{
  this->make_copy(another);
}

void MetropolisStandardMCMCOutput::operator=(const MetropolisStandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* MetropolisStandardMCMCOutput::duplicate() const
{
  return( new MetropolisStandardMCMCOutput(*this));
}

void MetropolisStandardMCMCOutput::make_copy(const MetropolisStandardMCMCOutput &another)
{
  if (another.mcmc!=NULL)
    this->mcmc = another.mcmc->metropolis_mcmc_duplicate();
  else
    this->mcmc = NULL;
}

MCMC* MetropolisStandardMCMCOutput::get_mcmc()
{
  return this->mcmc;
}

const MCMC* MetropolisStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
