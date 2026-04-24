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
  /**
   * @file truncated_gaussian_distribution_factor.h
   * @brief Defines the TruncatedGaussianDistributionFactor class.
   *
   * A truncated gaussian distribution factor. Evaluates the log-density of a truncated gaussian prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class TruncatedGaussianDistributionFactor
   * @brief A truncated gaussian distribution factor derived from DistributionFactor.
   */


class TruncatedGaussianDistributionFactor : public DistributionFactor
{
  
public:
  
  /**
   * @brief Default constructor for TruncatedGaussianDistributionFactor.
   */
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
  
  /**
   * @brief Destructor for TruncatedGaussianDistributionFactor.
   */
  virtual ~TruncatedGaussianDistributionFactor();
  
  /**
   * @brief Copy constructor for TruncatedGaussianDistributionFactor.
   *
   * @param another The TruncatedGaussianDistributionFactor instance to copy from.
   */
  TruncatedGaussianDistributionFactor(const TruncatedGaussianDistributionFactor &another);
  
  /**
   * @brief Assignment operator for TruncatedGaussianDistributionFactor.
   *
   * @param another The TruncatedGaussianDistributionFactor instance to copy from.
   */
  void operator=(const TruncatedGaussianDistributionFactor &another);
  /**
   * @brief Creates a deep copy of this TruncatedGaussianDistributionFactor object.
   *
   * @return The result.
   */
  Factor* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a distribution_factor pointer.
   *
   * @return The result.
   */
  DistributionFactor* distribution_factor_duplicate() const;
  
  /**
   * @brief Performs the distribution evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;
  
protected:
  
  /**
   * @brief Copies the state of another TruncatedGaussianDistributionFactor into this object.
   *
   * @param another The TruncatedGaussianDistributionFactor instance to copy from.
   */
  void make_copy(const TruncatedGaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
  boost::unordered_map< std::string, arma::colvec> lower_info;
  boost::unordered_map< std::string, arma::colvec> upper_info;
  
};
}

#endif
