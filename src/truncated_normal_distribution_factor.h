#ifndef TRUNCATEDNORMALDISTRIBUTIONFACTOR_H
#define TRUNCATEDNORMALDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "distributions.h" // (log_normal_mass, C218)

namespace ilike
{

// log density of the normal(mean, sd) truncated to [lower, upper] at x
inline double log_truncated_normal(double x, double mean, double sd, double lower, double upper)
{
  if ((x < lower) || (x > upper) || !(sd > 0.0) || !(lower < upper))
    return -arma::datum::inf;
  double z = (x - mean)/sd;
  return -0.5*z*z - 0.5*std::log(2.0*arma::datum::pi) - std::log(sd) -
         log_normal_mass((lower - mean)/sd, (upper - mean)/sd);
}

// The prior tnorm(x, mean, sd, lower, upper) (C217): the normal(mean, sd)
// truncated to [lower, upper] (either may be infinite), for one value or
// each element of a block independently.
class TruncatedNormalDistributionFactor : public DistributionFactor
{
public:

  TruncatedNormalDistributionFactor();
  TruncatedNormalDistributionFactor(const std::string &variable_in, double mean_in, double sd_in,
                                    double lower_in, double upper_in);
  virtual ~TruncatedNormalDistributionFactor();
  TruncatedNormalDistributionFactor(const TruncatedNormalDistributionFactor &another);
  void operator=(const TruncatedNormalDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;

  double distribution_evaluate(const Parameters &input) const;
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;

protected:

  void make_copy(const TruncatedNormalDistributionFactor &another);
  std::string variable;
  double mean;
  double sd;
  double lower;
  double upper;
};

}

#endif
