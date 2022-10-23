#ifndef GAUSSIANMCMCADAPTOR_H
#define GAUSSIANMCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "mcmc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"
#include "gaussian_proposal_info.h"

class SMCOutput;
class ScaleRecursiveParameterEstimator;
class VectorRecursiveParameterEstimator;
class GaussianRecursiveParameterEstimator;

class GaussianMCMCAdaptor : public MCMCAdaptor
{

public:

  GaussianMCMCAdaptor();
  virtual ~GaussianMCMCAdaptor();

  GaussianMCMCAdaptor(const GaussianMCMCAdaptor &another);

  void operator=(const GaussianMCMCAdaptor &another);
  MCMCAdaptor* duplicate() const;

protected:

  void make_copy(const GaussianMCMCAdaptor &another);
  
  void specific_mcmc_adapt(const Particle &latest_particle,
                           size_t iteration_counter);
  
  GaussianProposalInfo initial_proposal_info;
  
  std::vector<std::string> variable_names;
  
  // not stored here
  std::vector<GaussianProposalInfo*> gaussian_info_pointers;
  
  // stored here
  std::vector<ScaleRecursiveParameterEstimator*> scale_estimators;
  std::vector<VectorRecursiveParameterEstimator*> mean_estimators;
  std::vector<GaussianRecursiveParameterEstimator*> gaussian_estimators;

};

#endif
