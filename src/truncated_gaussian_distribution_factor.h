#ifndef TRUNCATEDGAUSSIANDISTRIBUTIONFACTOR_H
#define TRUNCATEDGAUSSIANDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
class TruncatedGaussianDistributionFactor : public DistributionFactor
{
  
public:
  
  TruncatedGaussianDistributionFactor();
  
  TruncatedGaussianDistributionFactor(const std::string &variable_in,
                                      double mean_in,
                                      double sd_in,
                                      double lower_in,
                                      double upper_in);
  
  TruncatedGaussianDistributionFactor(const std::string &variable_in,
                                      const arma::colvec &mean_in,
                                      const arma::mat &covariance_in,
                                      const arma::colvec &lower_in,
                                      const arma::colvec &upper_in);
  
  TruncatedGaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
                                      const std::vector<arma::colvec> &means_in,
                                      const std::vector<arma::mat> &covariances_in,
                                      const std::vector<arma::colvec> &lowers_in,
                                      const std::vector<arma::colvec> &uppers_in);
  
  virtual ~TruncatedGaussianDistributionFactor();
  
  TruncatedGaussianDistributionFactor(const TruncatedGaussianDistributionFactor &another);
  
  void operator=(const TruncatedGaussianDistributionFactor &another);
  Factor* duplicate() const;
  DistributionFactor* distribution_factor_duplicate() const;
  
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;
  
protected:
  
  void make_copy(const TruncatedGaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
  boost::unordered_map< std::string, arma::colvec> lower_info;
  boost::unordered_map< std::string, arma::colvec> upper_info;
  
};
}

#endif
