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
  // (C269, PMMH P1) the covariance is a constant (the language knows it): its
  // Cholesky factor, inverse and log-determinant are computed here, once,
  // instead of at every call (computed eagerly: no lazy first call that
  // parallel filter chunks could race on)
  void set_covariance_constant();
  // (C269, PMMH P1, owner's "declared reads") the matrix and covariance
  // functions read only these variables (the user's reads = [...]): at the
  // start of a filter run, prepare(parameters) evaluates them once on a
  // Parameters holding only those variables (so a function that reads any
  // other variable stops with a message: the declaration is checked) and
  // keeps their values with the variables' values; a later call uses the kept
  // values only if the particle's values of the declared variables are equal
  // to the kept ones, and otherwise computes as before (so a stale or missing
  // preparation costs time, never correctness). prepare is called at
  // sequential points only; the calls during the run only read.
  void set_declared_reads(const std::vector<std::string> &reads_in);
  void prepare(const Parameters &parameters);
  
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
  bool covariance_constant = false;
  std::vector<std::string> declared_reads;
  bool cache_valid = false;
  std::vector<arma::mat> cached_values; // (of the declared variables)
  arma::mat cached_A;
  bool cached_A_set = false;
  arma::mat cached_chol;
  arma::mat cached_inv;
  double cached_logdet = 0.0;
  bool cached_covariance_set = false;
  bool cache_matches(const Parameters &parameters) const;
  arma::mat constant_chol;
  arma::mat constant_inv;
  double constant_logdet = 0.0;
  
};
}

#endif
