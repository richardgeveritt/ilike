#ifndef DENSITYESTIMATOR_H
#define DENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

namespace ilike
{
//class SMCOutput;
class DensityEstimatorOutput;

class DensityEstimator
{
  
public:
  
  DensityEstimator();
  DensityEstimator(const std::vector<std::string> &variables_in);
  virtual ~DensityEstimator();
  
  DensityEstimator(const DensityEstimator &another);
  
  void operator=(const DensityEstimator &another);
  virtual DensityEstimator* duplicate() const=0;
  
  virtual DensityEstimatorOutput* initialise()=0;
  
  std::vector<std::string> get_variables() const;
  
protected:
  
  std::vector<std::string> variables;
  
  void make_copy(const DensityEstimator &another);
  
};
}

#endif
