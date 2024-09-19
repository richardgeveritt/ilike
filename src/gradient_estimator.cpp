#include "gradient_estimator.h"
#include "gradient_estimator_output.h"
#include "proposal_kernel.h"
#include "direct_gradient_estimator.h"

namespace ilike
{
GradientEstimator::GradientEstimator()
{
  this->proposal = NULL;
}

GradientEstimator::~GradientEstimator()
{
  
}

GradientEstimator::GradientEstimator(const ProposalKernel* proposal_in)
{
  this->proposal = proposal_in;
}

GradientEstimator::GradientEstimator(const GradientEstimator &another)
{
  this->make_copy(another);
}

void GradientEstimator::operator=(const GradientEstimator &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void GradientEstimator::make_copy(const GradientEstimator &another)
{
  this->proposal = another.proposal;
}

void GradientEstimator::set_proposal(const ProposalKernel* proposal_in)
{
  this->proposal = proposal_in;
}
}
