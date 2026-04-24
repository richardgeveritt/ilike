#ifndef LOGGAUSSIANDISTRIBUTIONFACTOR_H
#define LOGGAUSSIANDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file loggaussian_distribution_factor.h
   * @brief Defines the LogGaussianDistributionFactor class.
   *
   * A log gaussian distribution factor. Evaluates the log-density of a log gaussian prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class LogGaussianDistributionFactor
   * @brief A log gaussian distribution factor derived from DistributionFactor.
   */


class LogGaussianDistributionFactor : public DistributionFactor
{
  
public:
  
  /**
   * @brief Default constructor for LogGaussianDistributionFactor.
   */
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
  
  /**
   * @brief Destructor for LogGaussianDistributionFactor.
   */
  virtual ~LogGaussianDistributionFactor();
  
  /**
   * @brief Copy constructor for LogGaussianDistributionFactor.
   *
   * @param another The LogGaussianDistributionFactor instance to copy from.
   */
  LogGaussianDistributionFactor(const LogGaussianDistributionFactor &another);
  
  /**
   * @brief Assignment operator for LogGaussianDistributionFactor.
   *
   * @param another The LogGaussianDistributionFactor instance to copy from.
   */
  void operator=(const LogGaussianDistributionFactor &another);
  /**
   * @brief Creates a deep copy of this LogGaussianDistributionFactor object.
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
   * @brief Copies the state of another LogGaussianDistributionFactor into this object.
   *
   * @param another The LogGaussianDistributionFactor instance to copy from.
   */
  void make_copy(const LogGaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
};

}

#endif
