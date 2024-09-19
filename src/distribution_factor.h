#ifndef DISTRIBUTIONFACTOR_H
#define DISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "factor.h"

namespace ilike
{
class DistributionFactor : public Factor
{
  
public:
  
  DistributionFactor();
  virtual ~DistributionFactor();
  
  DistributionFactor(const DistributionFactor &another);
  
  void operator=(const DistributionFactor &another);
  virtual DistributionFactor* distribution_factor_duplicate() const=0;
  
  double evaluate(const Parameters &input) const;
  virtual double distribution_evaluate(const Parameters &input) const=0;
  
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &input) const;
  virtual arma::mat distribution_evaluate_gradient(const std::string &variable,
                                                   const Parameters &input) const=0;
  
protected:
  
  void make_copy(const DistributionFactor &another);
  
};
}

#endif
