#ifndef GAUSSIANSMCADAPTOR_H
#define GAUSSIANSMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "smc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"
#include "gaussian_proposal_info.h"

class SMCOutput;
class VectorParameterEstimator;
class MatrixParameterEstimator;

class GaussianSMCAdaptor : public SMCAdaptor
{

public:

  GaussianSMCAdaptor();
  virtual ~GaussianSMCAdaptor();

  GaussianSMCAdaptor(const GaussianSMCAdaptor &another);

  void operator=(const GaussianSMCAdaptor &another);
  SMCAdaptor* duplicate() const;
  
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);

protected:

  void make_copy(const GaussianSMCAdaptor &another);
  
  std::vector<std::string> variable_names;
  
  // not stored here
  std::vector<GaussianProposalInfo*> gaussian_info_pointers;
  
  // stored here
  std::vector<Scale> scales;
  std::vector<VectorParameterEstimator*> mean_estimators;
  std::vector<MatrixParameterEstimator*> covariance_estimators;

};

#endif
