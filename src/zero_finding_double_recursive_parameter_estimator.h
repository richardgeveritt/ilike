#ifndef ZEROFINDINGDOUBLEPARAMETERESTIMATOR_H
#define ZEROFINDINGDOUBLEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "double_recursive_parameter_estimator.h"
#include "ilike_header.h"

namespace ilike
{
class ZeroFindingDoubleRecursiveParameterEstimator : public DoubleRecursiveParameterEstimator
{
  
public:
  
  ZeroFindingDoubleRecursiveParameterEstimator();
  ZeroFindingDoubleRecursiveParameterEstimator(double initial_value,
                                               double target_score_in);
  
  virtual ~ZeroFindingDoubleRecursiveParameterEstimator();
  
  ZeroFindingDoubleRecursiveParameterEstimator(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  
  void operator=(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  RecursiveParameterEstimator* duplicate() const;
  DoubleRecursiveParameterEstimator* double_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  GainPtr gain;
  
  double target_score;
  
  void make_copy(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  
};
}

#endif
