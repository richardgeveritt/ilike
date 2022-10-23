#include "mcmc_adaptor.h"
#include "proposal_kernel.h"

MCMCAdaptor::MCMCAdaptor()
{
}

MCMCAdaptor::~MCMCAdaptor()
{
}

MCMCAdaptor::MCMCAdaptor(const MCMCAdaptor &another)
{
  this->make_copy(another);
}

void MCMCAdaptor::operator=(const MCMCAdaptor &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void MCMCAdaptor::make_copy(const MCMCAdaptor &another)
{
  this->proposal = another.proposal;
}

void MCMCAdaptor::mcmc_adapt(Particle &latest_particle,
                             size_t iteration_counter)
{
  this->proposal->use_transform(latest_particle);
  this->specific_mcmc_adapt(latest_particle,
                            iteration_counter);
}
