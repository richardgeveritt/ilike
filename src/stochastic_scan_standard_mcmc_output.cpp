#include "stochastic_scan_standard_mcmc_output.h"
#include "utils.h"
#include "stochastic_scan_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput()
:StandardMCMCOutput()
{
}

StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput(StochasticScanMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<StochasticScanMCMC, OptionalDelete>(mcmc_in, OptionalDelete(true));
  this->termination->set_parameters(this);
}

StochasticScanStandardMCMCOutput::StochasticScanStandardMCMCOutput(StochasticScanMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<StochasticScanMCMC, OptionalDelete>(mcmc_in, OptionalDelete(owns_mcmc_in));
  this->termination->set_parameters(this);
}

StochasticScanStandardMCMCOutput::~StochasticScanStandardMCMCOutput()
{
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

  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* StochasticScanStandardMCMCOutput::duplicate() const
{
  return( new StochasticScanStandardMCMCOutput(*this));
}

void StochasticScanStandardMCMCOutput::make_copy(const StochasticScanStandardMCMCOutput &another)
{
  if (another.mcmc==NULL)
    this->mcmc = std::unique_ptr<StochasticScanMCMC, OptionalDelete>(nullptr, OptionalDelete(true));
  else if (another.mcmc.get_deleter().owns)
    this->mcmc = std::unique_ptr<StochasticScanMCMC, OptionalDelete>(another.mcmc->stochastic_scan_mcmc_duplicate(), OptionalDelete(true));
  else
    this->mcmc = std::unique_ptr<StochasticScanMCMC, OptionalDelete>(another.mcmc.get(), OptionalDelete(false));
}

MCMC* StochasticScanStandardMCMCOutput::get_mcmc()
{
  return this->mcmc.get();
}

const MCMC* StochasticScanStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc.get();
}
}
