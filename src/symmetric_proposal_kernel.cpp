#include "symmetric_proposal_kernel.h"
#include "mcmc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "smc_adaptor.h"

SymmetricProposalKernel::SymmetricProposalKernel()
  :ProposalKernel()
{
}

SymmetricProposalKernel::~SymmetricProposalKernel()
{
}

SymmetricProposalKernel::SymmetricProposalKernel(const SymmetricProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void SymmetricProposalKernel::operator=(const SymmetricProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

void SymmetricProposalKernel::make_copy(const SymmetricProposalKernel &another)
{
}
