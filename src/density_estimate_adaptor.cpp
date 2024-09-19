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
}

void DensityEstimateAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
}
}
