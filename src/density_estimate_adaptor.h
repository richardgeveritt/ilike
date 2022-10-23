#ifndef DENSITYESTIMATEADAPTOR_H
#define DENSITYESTIMATEADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "smc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"

class SMCOutput;

class DensityEstimateAdaptor : public SMCAdaptor
{

public:

  DensityEstimateAdaptor();
  virtual ~DensityEstimateAdaptor();

  DensityEstimateAdaptor(const DensityEstimateAdaptor &another);

  void operator=(const DensityEstimateAdaptor &another);
  SMCAdaptor* duplicate() const;
  
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);

protected:

  void make_copy(const DensityEstimateAdaptor &another);

  // needs density estimator hierarchy
  
};

#endif
