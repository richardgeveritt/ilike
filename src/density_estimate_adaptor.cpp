#include "density_estimate_adaptor.h"
#include "utils.h"
#include "smc_output.h"

namespace ilike
{
DensityEstimateAdaptor::DensityEstimateAdaptor()
:SMCAdaptor()
{
}

DensityEstimateAdaptor::~DensityEstimateAdaptor()
{
}

DensityEstimateAdaptor::DensityEstimateAdaptor(const DensityEstimateAdaptor &another)
:SMCAdaptor(another)
{
  // copy
  this->make_copy(another);
}

void DensityEstimateAdaptor::operator=(const DensityEstimateAdaptor &another)
{
  if(this == &another)
    return;

  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

void DensityEstimateAdaptor::make_copy(const DensityEstimateAdaptor &another)
{
}

SMCAdaptor* DensityEstimateAdaptor::duplicate() const
{
  return( new DensityEstimateAdaptor(*this));
}

void DensityEstimateAdaptor::smc_adapt(SMCOutput* current_state)
{
  // reads the generation's particles one by one: not available for a
  // generation stored flat (dev/plan.md 2.6a; no adaptor is set from R)
  if (current_state->back().is_flat())
    Rcpp::stop("DensityEstimateAdaptor::smc_adapt: not available for a generation stored flat.");
}

void DensityEstimateAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
}
}
