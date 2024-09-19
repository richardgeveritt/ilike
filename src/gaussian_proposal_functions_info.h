#ifndef GAUSSIANPROPOSALFUNCTIONSINFO_H
#define GAUSSIANPROPOSALFUNCTIONSINFO_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "scale.h"
#include "ilike_header.h"

namespace ilike
{
class GaussianProposalFunctionsInfo
{
  
public:
  
  GaussianProposalFunctionsInfo();
  
  virtual ~GaussianProposalFunctionsInfo();
  
  GaussianProposalFunctionsInfo(const GaussianProposalFunctionsInfo &another);
  
  GaussianProposalFunctionsInfo& operator=(const GaussianProposalFunctionsInfo &another);
  GaussianProposalFunctionsInfo* duplicate() const;
  
  GaussianProposalFunctionsInfo(GaussianProposalFunctionsInfo &&another);
  GaussianProposalFunctionsInfo& operator=(GaussianProposalFunctionsInfo &&another);
  
  void set_mean(const GetVectorPtr &mean_in);
  void set_covariance(const GetMatrixPtr &covariance_in);
  void set_A(const GetMatrixPtr &A_in);
  //void set_only_covariance(const arma::mat &covariance_in);
  
  arma::colvec get_mean(const Parameters &parameters) const;
  arma::mat get_A(const Parameters &parameters) const;
  arma::mat get_covariance(const Parameters &parameters) const;
  double get_double_scale() const;
  Scale get_scale() const;
  
  arma::mat get_chol(const Parameters &parameters) const;
  arma::mat get_inv(const Parameters &parameters) const;
  double get_logdet(const Parameters &parameters) const;
  
  double& get_double_scale();
  Scale& get_scale();
  
protected:
  
  void make_copy(const GaussianProposalFunctionsInfo &another);
  
  void make_copy(GaussianProposalFunctionsInfo &&another);
  
  GetVectorPtr mean;
  GetMatrixPtr A;
  GetMatrixPtr covariance;
  Scale scale;
  
  arma::mat chol;
  
};
}

#endif
