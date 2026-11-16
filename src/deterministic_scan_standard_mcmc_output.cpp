#include "deterministic_scan_standard_mcmc_output.h"
#include "utils.h"
#include "deterministic_scan_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput()
:StandardMCMCOutput()
{
}

DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<DeterministicScanMCMC, OptionalDelete>(mcmc_in, OptionalDelete(true));
  this->termination->set_parameters(this);
}

DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = std::unique_ptr<DeterministicScanMCMC, OptionalDelete>(mcmc_in, OptionalDelete(owns_mcmc_in));
  this->termination->set_parameters(this);
}

DeterministicScanStandardMCMCOutput::~DeterministicScanStandardMCMCOutput()
{
}

//Copy constructor for the DeterministicScanStandardMCMCOutput class.
DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput(const DeterministicScanStandardMCMCOutput &another)
:StandardMCMCOutput(another)
{
  this->make_copy(another);
}

void DeterministicScanStandardMCMCOutput::operator=(const DeterministicScanStandardMCMCOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* DeterministicScanStandardMCMCOutput::duplicate() const
{
  return( new DeterministicScanStandardMCMCOutput(*this));
}

void DeterministicScanStandardMCMCOutput::make_copy(const DeterministicScanStandardMCMCOutput &another)
{
  if (another.mcmc==NULL)
    this->mcmc = std::unique_ptr<DeterministicScanMCMC, OptionalDelete>(nullptr, OptionalDelete(true));
  else if (another.mcmc.get_deleter().owns)
    this->mcmc = std::unique_ptr<DeterministicScanMCMC, OptionalDelete>(another.mcmc->deterministic_scan_mcmc_duplicate(), OptionalDelete(true));
  else
    this->mcmc = std::unique_ptr<DeterministicScanMCMC, OptionalDelete>(another.mcmc.get(), OptionalDelete(false));
}

/*
 void DeterministicScanStandardMCMCOutput::specific_mcmc_adapt()
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
 Particle DeterministicScanStandardMCMCOutput::move(RandomNumberGenerator &rng,
 const Particle &particle) const
 {
 Particle current_particle = particle;
 for (std::vector<size_t>::const_iterator i=this->order.begin();
 i!=this->order.end();
 ++i)
 {
 current_particle = this->moves[*i]->move(rng,
 current_particle);
 }
 return current_particle;
 }
 
 Particle DeterministicScanStandardMCMCOutput::subsample_move(RandomNumberGenerator &rng,
 const Particle &particle) const
 {
 Particle current_particle = particle;
 for (std::vector<size_t>::const_iterator i=this->order.begin();
 i!=this->order.end();
 ++i)
 {
 current_particle = this->moves[*i]->subsample_move(rng,
 current_particle);
 }
 return current_particle;
 }
 */

MCMC* DeterministicScanStandardMCMCOutput::get_mcmc()
{
  return this->mcmc.get();
}

const MCMC* DeterministicScanStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc.get();
}
}
