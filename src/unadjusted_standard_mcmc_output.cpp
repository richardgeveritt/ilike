#include "unadjusted_standard_mcmc_output.h"
#include "utils.h"
#include "unadjusted_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput()
:StandardMCMCOutput()
{
}

UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<UnadjustedMCMC, OptionalDelete>(mcmc_in, OptionalDelete(true));
  this->termination->set_parameters(this);
}

UnadjustedStandardMCMCOutput::UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<UnadjustedMCMC, OptionalDelete>(mcmc_in, OptionalDelete(owns_mcmc_in));
  this->termination->set_parameters(this);
}

UnadjustedStandardMCMCOutput::~UnadjustedStandardMCMCOutput()
{
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
  
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* UnadjustedStandardMCMCOutput::duplicate() const
{
  return( new UnadjustedStandardMCMCOutput(*this));
}

void UnadjustedStandardMCMCOutput::make_copy(const UnadjustedStandardMCMCOutput &another)
{
  if (another.mcmc==NULL)
    this->mcmc = std::unique_ptr<UnadjustedMCMC, OptionalDelete>(nullptr, OptionalDelete(true));
  else if (another.mcmc.get_deleter().owns)
    this->mcmc = std::unique_ptr<UnadjustedMCMC, OptionalDelete>(another.mcmc->unadjusted_mcmc_duplicate(), OptionalDelete(true));
  else
    this->mcmc = std::unique_ptr<UnadjustedMCMC, OptionalDelete>(another.mcmc.get(), OptionalDelete(false));
}

MCMC* UnadjustedStandardMCMCOutput::get_mcmc()
{
  return this->mcmc.get();
}

const MCMC* UnadjustedStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc.get();
}
}
