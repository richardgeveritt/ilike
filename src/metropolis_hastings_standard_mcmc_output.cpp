#include "metropolis_hastings_standard_mcmc_output.h"
#include "utils.h"
#include "metropolis_hastings_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput()
:StandardMCMCOutput()
{
}

MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<MetropolisHastingsMCMC, OptionalDelete>(mcmc_in, OptionalDelete(true));
  this->termination->set_parameters(this);
}

MetropolisHastingsStandardMCMCOutput::MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<MetropolisHastingsMCMC, OptionalDelete>(mcmc_in, OptionalDelete(owns_mcmc_in));
  this->termination->set_parameters(this);
}

MetropolisHastingsStandardMCMCOutput::~MetropolisHastingsStandardMCMCOutput()
{
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

  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* MetropolisHastingsStandardMCMCOutput::duplicate() const
{
  return( new MetropolisHastingsStandardMCMCOutput(*this));
}

void MetropolisHastingsStandardMCMCOutput::make_copy(const MetropolisHastingsStandardMCMCOutput &another)
{
  if (another.mcmc==NULL)
    this->mcmc = std::unique_ptr<MetropolisHastingsMCMC, OptionalDelete>(nullptr, OptionalDelete(true));
  else if (another.mcmc.get_deleter().owns)
    this->mcmc = std::unique_ptr<MetropolisHastingsMCMC, OptionalDelete>(another.mcmc->metropolis_hastings_mcmc_duplicate(), OptionalDelete(true));
  else
    this->mcmc = std::unique_ptr<MetropolisHastingsMCMC, OptionalDelete>(another.mcmc.get(), OptionalDelete(false));
}

MCMC* MetropolisHastingsStandardMCMCOutput::get_mcmc()
{
  return this->mcmc.get();
}

const MCMC* MetropolisHastingsStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc.get();
}
}
