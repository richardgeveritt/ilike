#ifndef GAUSSIANPROPOSALINFO_H
#define GAUSSIANPROPOSALINFO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "scale.h"

class GaussianProposalInfo
{

public:

  GaussianProposalInfo();

  virtual ~GaussianProposalInfo();

  GaussianProposalInfo(const arma::colvec &mean_in);
  GaussianProposalInfo(const arma::mat &covariance_in);
  GaussianProposalInfo(const arma::colvec &mean_in,
                       const arma::mat &covariance_in);
  
  GaussianProposalInfo(const GaussianProposalInfo &another);

  GaussianProposalInfo& operator=(const GaussianProposalInfo &another);
  GaussianProposalInfo* duplicate() const;
  
  GaussianProposalInfo(GaussianProposalInfo &&another);
  GaussianProposalInfo& operator=(GaussianProposalInfo &&another);
  
  void set_mean(const arma::colvec &mean_in);
  void set_covariance(const arma::mat &covariance_in);
  //void set_only_covariance(const arma::mat &covariance_in);
  void set_covariance_info();
  
  arma::colvec get_mean() const;
  arma::mat get_covariance() const;
  double get_double_scale() const;
  Scale get_scale() const;
  arma::mat get_chol() const;
  arma::mat get_inv() const;
  arma::mat get_inv_chol() const;
  double get_logdet() const;
  
  arma::colvec& get_mean();
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
  arma::mat covariance;
  Scale scale;
  arma::mat chol;
  arma::mat inv;
  arma::mat inv_chol;
  double logdet;

};

#endif
