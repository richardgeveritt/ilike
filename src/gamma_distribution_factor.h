#ifndef GAMMADISTRIBUTIONFACTOR_H
#define GAMMADISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"

class GammaDistributionFactor : public DistributionFactor
{

public:

  GammaDistributionFactor();

  GammaDistributionFactor(const std::string &variable_in,
                          double shape_in,
                          double rate_in);
  
  virtual ~GammaDistributionFactor();

  GammaDistributionFactor(const GammaDistributionFactor &another);

  void operator=(const GammaDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;
  
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;

protected:

  void make_copy(const GammaDistributionFactor &another);
  
  std::string variable;
  double shape;
  double rate;

};

#endif
