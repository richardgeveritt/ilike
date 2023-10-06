#ifndef GAUSSIANNOISEPROPOSALKERNEL_H
#define GAUSSIANNOISEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

class GaussianNoiseProposalKernel : public ProposalKernel
{

public:

  GaussianNoiseProposalKernel();
  virtual ~GaussianNoiseProposalKernel();

  GaussianNoiseProposalKernel(const GaussianNoiseProposalKernel &another);

  void operator=(const GaussianNoiseProposalKernel &another);
  //Kernel* duplicate() const;
  virtual GaussianNoiseProposalKernel* gaussian_noise_proposal_kernel_duplicate() const=0;

protected:

  void make_copy(const GaussianNoiseProposalKernel &another);
  
};

#endif
