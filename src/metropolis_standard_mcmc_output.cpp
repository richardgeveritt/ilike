#include "metropolis_standard_mcmc_output.h"
#include "utils.h"
#include "metropolis_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
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

MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->owns_mcmc = owns_mcmc_in;
  this->termination->set_parameters(this);
}

MetropolisStandardMCMCOutput::~MetropolisStandardMCMCOutput()
{
  if ((this->mcmc!=NULL) && this->owns_mcmc)
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
  
  if ((this->mcmc!=NULL) && this->owns_mcmc)
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
  this->owns_mcmc = another.owns_mcmc;
  if (another.mcmc==NULL)
    this->mcmc = NULL;
  else if (another.owns_mcmc)
    this->mcmc = another.mcmc->metropolis_mcmc_duplicate();
  else
    this->mcmc = another.mcmc;
}

MCMC* MetropolisStandardMCMCOutput::get_mcmc()
{
  return this->mcmc;
}

const MCMC* MetropolisStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
}
