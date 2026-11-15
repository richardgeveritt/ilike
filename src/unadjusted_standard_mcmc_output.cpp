#include "unadjusted_standard_mcmc_output.h"
#include "utils.h"
#include "unadjusted_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput()
:StandardMCMCOutput()
{
  this->mcmc = NULL;
}

UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->termination->set_parameters(this);
}

UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->owns_mcmc = owns_mcmc_in;
  this->termination->set_parameters(this);
}

UnadjustedStandardMCMCOutput::~UnadjustedStandardMCMCOutput()
{
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
}

//Copy constructor for the UnadjustedStandardMCMCOutput class.
UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput(const UnadjustedStandardMCMCOutput &another)
:StandardMCMCOutput(another)
{
  this->make_copy(another);
}

void UnadjustedStandardMCMCOutput::operator=(const UnadjustedStandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* UnadjustedStandardMCMCOutput::duplicate() const
{
  return( new UnadjustedStandardMCMCOutput(*this));
}

void UnadjustedStandardMCMCOutput::make_copy(const UnadjustedStandardMCMCOutput &another)
{
  this->owns_mcmc = another.owns_mcmc;
  if (another.mcmc==NULL)
    this->mcmc = NULL;
  else if (another.owns_mcmc)
    this->mcmc = another.mcmc->unadjusted_mcmc_duplicate();
  else
    this->mcmc = another.mcmc;
}

MCMC* UnadjustedStandardMCMCOutput::get_mcmc()
{
  return this->mcmc;
}

const MCMC* UnadjustedStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
}
