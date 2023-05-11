#ifndef LOGGAUSSIANDISTRIBUTIONFACTOR_H
#define LOGGAUSSIANDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class LogGaussianDistributionFactor : public DistributionFactor
{
  
public:
  
  LogGaussianDistributionFactor();
  
  LogGaussianDistributionFactor(const std::string &variable_name_in,
                                const double &mean_in,
                                const double &covariance_in);
  
  LogGaussianDistributionFactor(const std::string &variable_name_in,
                                const arma::colvec &mean_in,
                                const arma::mat &covariance_in);
  
  LogGaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
                                const std::vector<arma::colvec> &means_in,
                                const std::vector<arma::mat> &covariances_in);
  
  virtual ~LogGaussianDistributionFactor();
  
  LogGaussianDistributionFactor(const LogGaussianDistributionFactor &another);
  
  void operator=(const LogGaussianDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;
  
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;
  
protected:
  
  void make_copy(const LogGaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
};

#endif
