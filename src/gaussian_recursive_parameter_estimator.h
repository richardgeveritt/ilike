#ifndef GAUSSIANRECURSIVEPARAMETERESTIMATOR_H
#define GAUSSIANRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include <boost/unordered_map.hpp>
#include "recursive_parameter_estimator.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
class GaussianRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  GaussianRecursiveParameterEstimator();
  
  virtual ~GaussianRecursiveParameterEstimator();
  
  GaussianRecursiveParameterEstimator(const GaussianRecursiveParameterEstimator &another);
  
  void operator=(const GaussianRecursiveParameterEstimator &another);
  virtual GaussianRecursiveParameterEstimator* gaussian_duplicate() const=0;
  
  GaussianProposalInfo estimated;
  
protected:
  
  void make_copy(const GaussianRecursiveParameterEstimator &another);
  
};
}

#endif
