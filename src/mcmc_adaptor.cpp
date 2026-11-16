#include "mcmc_adaptor.h"
#include "proposal_kernel.h"

namespace ilike
{

MCMCAdaptor::MCMCAdaptor()
{
  this->proposal = NULL;
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
  // the copy belongs to a copy of the kernel, which points it at itself
  // (ProposalKernel::make_copy); until then it points at the original's
  this->proposal = another.proposal;
}

void MCMCAdaptor::set_proposal(ProposalKernel* proposal_in)
{
  this->proposal = proposal_in;
}

void MCMCAdaptor::prepare_for_move()
{
}

void MCMCAdaptor::mcmc_adapt(const Particle &latest_particle,
                             size_t iteration_counter)
{
  Rcpp::stop("MCMC adaptation is available for Metropolis (m_proposal) and Metropolis-Hastings (mh_proposal) moves only.");
}

}
