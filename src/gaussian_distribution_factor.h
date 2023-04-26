#ifndef GAUSSIANDISTRIBUTIONFACTOR_H
#define GAUSSIANDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class GaussianDistributionFactor : public DistributionFactor
{

public:

  GaussianDistributionFactor();

  GaussianDistributionFactor(const std::string &variable_in,
                             double mean_in,
                             double sd_in);
  
  GaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
                             const std::vector<arma::colvec> &means_in,
                             const std::vector<arma::mat> &covariances_in);
  
  virtual ~GaussianDistributionFactor();

  GaussianDistributionFactor(const GaussianDistributionFactor &another);

  void operator=(const GaussianDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;
  
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;

protected:

  void make_copy(const GaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;

};

#endif
