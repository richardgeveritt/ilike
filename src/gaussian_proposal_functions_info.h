#ifndef GAUSSIANPROPOSALFUNCTIONSINFO_H
#define GAUSSIANPROPOSALFUNCTIONSINFO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "scale.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file gaussian_proposal_functions_info.h
   * @brief Defines the GaussianProposalFunctionsInfo class.
   *
   * Provides gaussian proposal functions info functionality.
   *
   * @namespace ilike
   * @class GaussianProposalFunctionsInfo
   * @brief The gaussian proposal functions info class.
   */


class GaussianProposalFunctionsInfo
{
  
public:
  
  /**
   * @brief Default constructor for GaussianProposalFunctionsInfo.
   */
  GaussianProposalFunctionsInfo();
  
  /**
   * @brief Destructor for GaussianProposalFunctionsInfo.
   */
  virtual ~GaussianProposalFunctionsInfo();
  
  /**
   * @brief Copy constructor for GaussianProposalFunctionsInfo.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  GaussianProposalFunctionsInfo(const GaussianProposalFunctionsInfo &another);
  
  /**
   * @brief Assignment operator for GaussianProposalFunctionsInfo.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  GaussianProposalFunctionsInfo& operator=(const GaussianProposalFunctionsInfo &another);
  /**
   * @brief Creates a deep copy of this GaussianProposalFunctionsInfo object.
   *
   * @return The result.
   */
  GaussianProposalFunctionsInfo* duplicate() const;
  
  /**
   * @brief Copy constructor for GaussianProposalFunctionsInfo.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  GaussianProposalFunctionsInfo(GaussianProposalFunctionsInfo &&another);
  /**
   * @brief Assignment operator for GaussianProposalFunctionsInfo.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  GaussianProposalFunctionsInfo& operator=(GaussianProposalFunctionsInfo &&another);
  
  /**
   * @brief Sets the mean.
   *
   * @param mean_in The mean.
   */
  void set_mean(const GetVectorPtr &mean_in);
  /**
   * @brief Sets the covariance.
   *
   * @param covariance_in The covariance.
   */
  void set_covariance(const GetMatrixPtr &covariance_in);
  /**
   * @brief Sets the a.
   *
   * @param A_in The a.
   */
  void set_A(const GetMatrixPtr &A_in);
  //void set_only_covariance(const arma::mat &covariance_in);
  
  /**
   * @brief Returns the mean.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::colvec get_mean(const Parameters &parameters) const;
  /**
   * @brief Returns the a.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::mat get_A(const Parameters &parameters) const;
  /**
   * @brief Returns the covariance.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::mat get_covariance(const Parameters &parameters) const;
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
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::mat get_chol(const Parameters &parameters) const;
  /**
   * @brief Returns the inv.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::mat get_inv(const Parameters &parameters) const;
  /**
   * @brief Returns the logdet.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double get_logdet(const Parameters &parameters) const;
  
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
  
protected:
  
  /**
   * @brief Copies the state of another GaussianProposalFunctionsInfo into this object.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  void make_copy(const GaussianProposalFunctionsInfo &another);
  
  /**
   * @brief Copies the state of another GaussianProposalFunctionsInfo into this object.
   *
   * @param another The GaussianProposalFunctionsInfo instance to copy from.
   */
  void make_copy(GaussianProposalFunctionsInfo &&another);
  
  /** @brief The mean. */
  GetVectorPtr mean;
  /** @brief The a. */
  GetMatrixPtr A;
  /** @brief The covariance. */
  GetMatrixPtr covariance;
  /** @brief The scale. */
  Scale scale;
  
  /** @brief The chol. */
  arma::mat chol;
  
};
}

#endif
