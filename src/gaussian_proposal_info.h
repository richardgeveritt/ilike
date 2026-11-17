#ifndef GAUSSIANPROPOSALINFO_H
#define GAUSSIANPROPOSALINFO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "scale.h"

namespace ilike
{

class GaussianProposalInfo
{
  
public:
  
  GaussianProposalInfo();
  
  virtual ~GaussianProposalInfo();
  
  GaussianProposalInfo(const double &mean_in,
                       const double &sd_in);
  
  GaussianProposalInfo(const double &sd_in);
  
  GaussianProposalInfo(const arma::colvec &mean_in);
  GaussianProposalInfo(const arma::mat &covariance_in);
  GaussianProposalInfo(const arma::mat &covariance_in,
                       double scale_in);
  GaussianProposalInfo(const arma::colvec &mean_in,
                       const arma::mat &covariance_in);
  GaussianProposalInfo(const arma::colvec &mean_in,
                       const arma::mat &covariance_in,
                       double scale_in);
  
  GaussianProposalInfo(const GaussianProposalInfo &another);
  
  GaussianProposalInfo& operator=(const GaussianProposalInfo &another);
  GaussianProposalInfo* duplicate() const;
  
  GaussianProposalInfo(GaussianProposalInfo &&another);
  GaussianProposalInfo& operator=(GaussianProposalInfo &&another);
  
  void set_mean(const arma::colvec &mean_in);
  void set_covariance(const arma::mat &covariance_in);
  void set_A(const arma::mat &A_in);
  void set_covariance_info();
  
  const arma::colvec& get_mean() const;
  arma::mat get_A() const;
  const arma::mat& get_covariance() const;
  double get_double_scale() const;
  Scale get_scale() const;
  const arma::mat& get_chol() const;
  // Lower-triangular factor L = chol^T (covariance = L L^T), kept so that
  // simulation can use L z: with R's reference BLAS a product with a
  // transposed matrix (chol^T z) is about 3 times slower at d = 1000.
  const arma::mat& get_chol_lower() const;
  const arma::mat& get_inv() const;
  const arma::mat& get_inv_chol() const;
  double get_logdet() const;
  
  arma::colvec& get_mean();
  arma::mat& get_A();
  arma::mat& get_covariance();
  double& get_double_scale();
  Scale& get_scale();
  arma::mat& get_chol();
  arma::mat& get_inv();
  arma::mat& get_inv_chol();
  double& get_logdet();
  
protected:
  
  void make_copy(const GaussianProposalInfo &another);
  
  void make_copy(GaussianProposalInfo &&another);
  
  arma::colvec mean;
  arma::mat A;
  arma::mat covariance;
  Scale scale;
  arma::mat chol;

  // chol^T, see get_chol_lower.
  arma::mat chol_lower;
  arma::mat inv;
  arma::mat inv_chol;
  double logdet;
  
};
}

#endif
