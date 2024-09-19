#ifndef CUSTOMDISTRIBUTIONFACTOR_H
#define CUSTOMDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"

namespace ilike
{
class CustomDistributionFactor : public DistributionFactor
{
  
public:
  
  CustomDistributionFactor();
  
  CustomDistributionFactor(EvaluateLogDistributionPtr distribution_in);
  
  CustomDistributionFactor(EvaluateLogDistributionPtr distribution_in,
                           EvaluateGradientLogDistributionPtr distribution_gradient_in);
  
  virtual ~CustomDistributionFactor();
  
  CustomDistributionFactor(const CustomDistributionFactor &another);
  
  void operator=(const CustomDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;
  
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;
  
protected:
  
  void make_copy(const CustomDistributionFactor &another);
  
  EvaluateLogDistributionPtr distribution;
  EvaluateGradientLogDistributionPtr distribution_gradient;
  
};
}

#endif
