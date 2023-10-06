#include "gaussian_noise_proposal_kernel.h"
#include "mcmc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "smc_adaptor.h"

GaussianNoiseProposalKernel::GaussianNoiseProposalKernel()
  :ProposalKernel()
{
}

GaussianNoiseProposalKernel::~GaussianNoiseProposalKernel()
{
}

GaussianNoiseProposalKernel::GaussianNoiseProposalKernel(const GaussianNoiseProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void GaussianNoiseProposalKernel::operator=(const GaussianNoiseProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

void GaussianNoiseProposalKernel::make_copy(const GaussianNoiseProposalKernel &another)
{
}
