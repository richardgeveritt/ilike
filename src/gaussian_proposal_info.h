#ifndef GAUSSIANPROPOSALINFO_H
#define GAUSSIANPROPOSALINFO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "scale.h"

namespace ilike
{
  /**
   * @file gaussian_proposal_info.h
   * @brief Defines the GaussianProposalInfo class.
   *
   * Provides gaussian proposal info functionality.
   *
   * @namespace ilike
   * @class GaussianProposalInfo
   * @brief The gaussian proposal info class.
   */


class GaussianProposalInfo
{
  
public:
  
  /**
   * @brief Default constructor for GaussianProposalInfo.
   */
  GaussianProposalInfo();
  
  /**
   * @brief Destructor for GaussianProposalInfo.
   */
  virtual ~GaussianProposalInfo();
  
  GaussianProposalInfo(const double &mean_in,
                       const double &sd_in);
  
  /**
   * @brief Constructs a GaussianProposalInfo object.
   *
   * @param sd_in The sd.
   */
  GaussianProposalInfo(const double &sd_in);
  
  /**
   * @brief Constructs a GaussianProposalInfo object.
   *
   * @param mean_in The mean.
   */
  GaussianProposalInfo(const arma::colvec &mean_in);
  /**
   * @brief Constructs a GaussianProposalInfo object.
   *
   * @param covariance_in The covariance.
   */
  GaussianProposalInfo(const arma::mat &covariance_in);
  GaussianProposalInfo(const arma::mat &covariance_in,
                       double scale_in);
  GaussianProposalInfo(const arma::colvec &mean_in,
                       const arma::mat &covariance_in);
  GaussianProposalInfo(const arma::colvec &mean_in,
                       const arma::mat &covariance_in,
                       double scale_in);
  
  /**
   * @brief Copy constructor for GaussianProposalInfo.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  GaussianProposalInfo(const GaussianProposalInfo &another);
  
  /**
   * @brief Assignment operator for GaussianProposalInfo.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  GaussianProposalInfo& operator=(const GaussianProposalInfo &another);
  /**
   * @brief Creates a deep copy of this GaussianProposalInfo object.
   *
   * @return The result.
   */
  GaussianProposalInfo* duplicate() const;
  
  /**
   * @brief Copy constructor for GaussianProposalInfo.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  GaussianProposalInfo(GaussianProposalInfo &&another);
  /**
   * @brief Assignment operator for GaussianProposalInfo.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  GaussianProposalInfo& operator=(GaussianProposalInfo &&another);
  
  /**
   * @brief Sets the mean.
   *
   * @param mean_in The mean.
   */
  void set_mean(const arma::colvec &mean_in);
  /**
   * @brief Sets the covariance.
   *
   * @param covariance_in The covariance.
   */
  void set_covariance(const arma::mat &covariance_in);
  /**
   * @brief Sets the a.
   *
   * @param A_in The a.
   */
  void set_A(const arma::mat &A_in);
  //void set_only_covariance(const arma::mat &covariance_in);
  /**
   * @brief Sets the covariance info.
   */
  void set_covariance_info();
  
  /**
   * @brief Returns the mean.
   *
   * @return The result.
   */
  arma::colvec get_mean() const;
  /**
   * @brief Returns the a.
   *
   * @return The result.
   */
  arma::mat get_A() const;
  /**
   * @brief Returns the covariance.
   *
   * @return The result.
   */
  arma::mat get_covariance() const;
  /**
   * @brief Returns the double scale.
   *
   * @return The result.
   */
  double get_double_scale() const;
  /**
   * @brief Returns the scale.
   *
   * @return The result.
   */
  Scale get_scale() const;
  /**
   * @brief Returns the chol.
   *
   * @return The result.
   */
  arma::mat get_chol() const;
  /**
   * @brief Returns the inv.
   *
   * @return The result.
   */
  arma::mat get_inv() const;
  /**
   * @brief Returns the inv chol.
   *
   * @return The result.
   */
  arma::mat get_inv_chol() const;
  /**
   * @brief Returns the logdet.
   *
   * @return The result.
   */
  double get_logdet() const;
  
  /**
   * @brief Returns the mean.
   *
   * @return The result.
   */
  arma::colvec& get_mean();
  /**
   * @brief Returns the a.
   *
   * @return The result.
   */
  arma::mat& get_A();
  /**
   * @brief Returns the covariance.
   *
   * @return The result.
   */
  arma::mat& get_covariance();
  /**
   * @brief Returns the double scale.
   *
   * @return The result.
   */
  double& get_double_scale();
  /**
   * @brief Returns the scale.
   *
   * @return The result.
   */
  Scale& get_scale();
  /**
   * @brief Returns the chol.
   *
   * @return The result.
   */
  arma::mat& get_chol();
  /**
   * @brief Returns the inv.
   *
   * @return The result.
   */
  arma::mat& get_inv();
  /**
   * @brief Returns the inv chol.
   *
   * @return The result.
   */
  arma::mat& get_inv_chol();
  /**
   * @brief Returns the logdet.
   *
   * @return The result.
   */
  double& get_logdet();
  
protected:
  
  /**
   * @brief Copies the state of another GaussianProposalInfo into this object.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  void make_copy(const GaussianProposalInfo &another);
  
  /**
   * @brief Copies the state of another GaussianProposalInfo into this object.
   *
   * @param another The GaussianProposalInfo instance to copy from.
   */
  void make_copy(GaussianProposalInfo &&another);
  
  /** @brief The mean. */
  arma::colvec mean;
  /** @brief The a. */
  arma::mat A;
  /** @brief The covariance. */
  arma::mat covariance;
  /** @brief The scale. */
  Scale scale;
  /** @brief The chol. */
  arma::mat chol;
  /** @brief The inv. */
  arma::mat inv;
  /** @brief The inv chol. */
  arma::mat inv_chol;
  /** @brief The logdet. */
  double logdet;
  
};
}

#endif
