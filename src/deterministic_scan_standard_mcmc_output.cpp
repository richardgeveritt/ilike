#include "deterministic_scan_standard_mcmc_output.h"
#include "utils.h"
#include "deterministic_scan_mcmc.h"
#include "mcmc_termination.h"

namespace ilike
{
DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput()
:StandardMCMCOutput()
{
  this->mcmc = NULL;
}

DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->termination->set_parameters(this);
}

DeterministicScanStandardMCMCOutput::DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in, bool owns_mcmc_in)
:StandardMCMCOutput(mcmc_in->get_duplicated_termination())
{
  this->mcmc = mcmc_in;
  this->owns_mcmc = owns_mcmc_in;
  this->termination->set_parameters(this);
}

DeterministicScanStandardMCMCOutput::~DeterministicScanStandardMCMCOutput()
{
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
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
  
  if ((this->mcmc!=NULL) && this->owns_mcmc)
    delete this->mcmc;
  
  StandardMCMCOutput::operator=(another);
  this->make_copy(another);
}

MoveOutput* DeterministicScanStandardMCMCOutput::duplicate() const
{
  return( new DeterministicScanStandardMCMCOutput(*this));
}

void DeterministicScanStandardMCMCOutput::make_copy(const DeterministicScanStandardMCMCOutput &another)
{
  this->owns_mcmc = another.owns_mcmc;
  if (another.mcmc==NULL)
    this->mcmc = NULL;
  else if (another.owns_mcmc)
    this->mcmc = another.mcmc->deterministic_scan_mcmc_duplicate();
  else
    this->mcmc = another.mcmc;
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
  return this->mcmc;
}

const MCMC* DeterministicScanStandardMCMCOutput::get_mcmc() const
{
  return this->mcmc;
}
}
