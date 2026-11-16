#include "metropolis_standard_mcmc_output.h"
#include "utils.h"
#include "metropolis_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput()
:StandardMCMCOutput()
{
}

MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<MetropolisMCMC, OptionalDelete>(mcmc_in, OptionalDelete(true));
  this->termination->set_parameters(this);
}

MetropolisStandardMCMCOutput::MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<MetropolisMCMC, OptionalDelete>(mcmc_in, OptionalDelete(owns_mcmc_in));
  this->termination->set_parameters(this);
}

MetropolisStandardMCMCOutput::~MetropolisStandardMCMCOutput()
{
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
  
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* MetropolisStandardMCMCOutput::duplicate() const
{
  return( new MetropolisStandardMCMCOutput(*this));
}

void MetropolisStandardMCMCOutput::make_copy(const MetropolisStandardMCMCOutput &another)
{
  if (another.mcmc==NULL)
    this->mcmc = std::unique_ptr<MetropolisMCMC, OptionalDelete>(nullptr, OptionalDelete(true));
  else if (another.mcmc.get_deleter().owns)
    this->mcmc = std::unique_ptr<MetropolisMCMC, OptionalDelete>(another.mcmc->metropolis_mcmc_duplicate(), OptionalDelete(true));
  else
    this->mcmc = std::unique_ptr<MetropolisMCMC, OptionalDelete>(another.mcmc.get(), OptionalDelete(false));
}

MCMC* MetropolisStandardMCMCOutput::get_mcmc()
{
  return this->mcmc.get();
}

const MCMC* MetropolisStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc.get();
}
}
