#ifndef GAUSSIANDISTRIBUTIONFACTOR_H
#define GAUSSIANDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file gaussian_distribution_factor.h
   * @brief Defines the GaussianDistributionFactor class.
   *
   * A gaussian distribution factor. Evaluates the log-density of a gaussian prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class GaussianDistributionFactor
   * @brief A gaussian distribution factor derived from DistributionFactor.
   */


class GaussianDistributionFactor : public DistributionFactor
{
  
public:
  
  /**
   * @brief Default constructor for GaussianDistributionFactor.
   */
  GaussianDistributionFactor();
  
  GaussianDistributionFactor(const std::string &variable_in,
                             double mean_in,
                             double sd_in);
  
  GaussianDistributionFactor(const std::string &variable_in,
                             const arma::colvec &mean_in,
                             const arma::mat &covariance_in);
  
  GaussianDistributionFactor(const std::vector<std::string> &variable_names_in,
                             const std::vector<arma::colvec> &means_in,
                             const std::vector<arma::mat> &covariances_in);
  
  /**
   * @brief Destructor for GaussianDistributionFactor.
   */
  virtual ~GaussianDistributionFactor();
  
  /**
   * @brief Copy constructor for GaussianDistributionFactor.
   *
   * @param another The GaussianDistributionFactor instance to copy from.
   */
  GaussianDistributionFactor(const GaussianDistributionFactor &another);
  
  /**
   * @brief Assignment operator for GaussianDistributionFactor.
   *
   * @param another The GaussianDistributionFactor instance to copy from.
   */
  void operator=(const GaussianDistributionFactor &another);
  /**
   * @brief Creates a deep copy of this GaussianDistributionFactor object.
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
   * @brief Copies the state of another GaussianDistributionFactor into this object.
   *
   * @param another The GaussianDistributionFactor instance to copy from.
   */
  void make_copy(const GaussianDistributionFactor &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
};
}

#endif
