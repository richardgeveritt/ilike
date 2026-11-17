#include "smc_adaptor.h"
#include "proposal_kernel.h"

namespace ilike
{

SMCAdaptor::SMCAdaptor()
{
  this->proposal = NULL;
}

SMCAdaptor::~SMCAdaptor()
{
}

SMCAdaptor::SMCAdaptor(const SMCAdaptor &another)
{
  this->make_copy(another);
}

void SMCAdaptor::operator=(const SMCAdaptor &another)
{
  if(this == &another)
    return;
  this->make_copy(another);
}

void SMCAdaptor::make_copy(const SMCAdaptor &another)
{
  // the copy belongs to a copy of the kernel, which points it at itself
  this->proposal = another.proposal;
}

void SMCAdaptor::set_proposal(ProposalKernel* proposal_in)
{
  this->proposal = proposal_in;
}

}
