#include "gaussian_proposal_info.h"
#include "utils.h"

GaussianProposalInfo::GaussianProposalInfo()
{
}

GaussianProposalInfo::GaussianProposalInfo(const double &mean_in,
                                           const double &sd_in)
: scale(1.0)
{
  this->mean = mean_in;
  arma::mat cov(1, 1, arma::fill::value(pow(sd_in,2.0)));
  this->set_covariance(cov);
}

GaussianProposalInfo::GaussianProposalInfo(const arma::mat &covariance_in)
: scale(1.0)
{
  this->set_covariance(covariance_in);
}

GaussianProposalInfo::GaussianProposalInfo(const arma::colvec &mean_in)
: scale(1.0)
{
  this->mean = mean_in;
}

GaussianProposalInfo::GaussianProposalInfo(const arma::colvec &mean_in,
                                           const arma::mat &covariance_in)
: scale(1.0)
{
  this->mean = mean_in;
  this->set_covariance(covariance_in);
}

GaussianProposalInfo::~GaussianProposalInfo()
{
  
}

//Copy constructor for the GaussianProposalInfo class.
GaussianProposalInfo::GaussianProposalInfo(const GaussianProposalInfo &another)
{
  this->make_copy(another);
}

GaussianProposalInfo& GaussianProposalInfo::operator=(const GaussianProposalInfo &another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(another);
  return *this;
}

GaussianProposalInfo* GaussianProposalInfo::duplicate() const
{
  return( new GaussianProposalInfo(*this));
}

void GaussianProposalInfo::make_copy(const GaussianProposalInfo &another)
{
  this->covariance = another.covariance;
  this->chol = another.chol;
  this->inv = another.inv;
  this->logdet = another.logdet;
  this->inv_chol = another.inv_chol;
  this->mean = another.mean;
  this->scale = another.scale;
}

//Move constructor for the GaussianProposalInfo class.
GaussianProposalInfo::GaussianProposalInfo(GaussianProposalInfo &&another)
: mean(), covariance(), scale(), chol(), inv(), inv_chol(), logdet()
{
  this->make_copy(std::move(another));
}

GaussianProposalInfo& GaussianProposalInfo::operator=(GaussianProposalInfo &&another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(std::move(another));
  return *this;
}

void GaussianProposalInfo::make_copy(GaussianProposalInfo &&another)
{
  this->covariance = std::move(another.covariance);
  this->chol = std::move(another.chol);
  this->inv = std::move(another.inv);
  this->logdet = std::move(another.logdet);
  this->inv_chol = std::move(another.inv_chol);
  this->mean = std::move(another.mean);
  this->scale = std::move(another.scale);
  
  another.covariance = arma::mat();
  another.chol = arma::mat();
  another.inv = arma::mat();
  another.logdet = 0.0;
  another.inv_chol = arma::mat();
  another.mean = arma::colvec();
  another.scale = Scale();
}

void GaussianProposalInfo::set_mean(const arma::colvec &mean_in)
{
  this->mean = mean_in;
}

void GaussianProposalInfo::set_covariance(const arma::mat &covariance_in)
{
  this->covariance = covariance_in;
  if ( (this->covariance.n_rows!=0) && (this->covariance.n_cols!=0) )
  {
    this->covariance = (this->covariance + this->covariance.t())/2.0;
    if (!this->covariance.is_sympd())
      Rcpp::stop("GaussianRandomWalkProposalKernel::setup_covariances - covariance matrix is not sympd.");
    this->chol= arma::chol(this->covariance);
    this->inv_chol = arma::inv(this->chol);
    this->inv = this->inv_chol*arma::trans(this->inv_chol);
    this->logdet = chol2logdet(this->chol);
  }
}

void GaussianProposalInfo::set_covariance_info()
{
  if ( (this->covariance.n_rows!=0) && (this->covariance.n_cols!=0) )
  {
    this->covariance = (this->covariance + this->covariance.t())/2.0;
    if (!this->covariance.is_sympd())
      Rcpp::stop("GaussianRandomWalkProposalKernel::setup_covariances - covariance matrix is not sympd.");
    this->chol= arma::chol(this->covariance);
    this->inv_chol = arma::inv(this->chol);
    this->inv = this->inv_chol*arma::trans(this->inv_chol);
    this->logdet = chol2logdet(this->chol);
  }
}

arma::colvec GaussianProposalInfo::get_mean() const
{
  return this->mean;
}

arma::mat GaussianProposalInfo::get_covariance() const
{
  return this->covariance;
}

double GaussianProposalInfo::get_double_scale() const
{
  return this->scale.get_constant();
}

Scale GaussianProposalInfo::get_scale() const
{
  return this->scale;
}

arma::mat GaussianProposalInfo::get_chol() const
{
  return this->chol;
}

arma::mat GaussianProposalInfo::get_inv() const
{
  return this->inv;
}

arma::mat GaussianProposalInfo::get_inv_chol() const
{
  return this->inv_chol;
}

double GaussianProposalInfo::get_logdet() const
{
  return this->logdet;
}

arma::colvec& GaussianProposalInfo::get_mean()
{
  return this->mean;
}

arma::mat& GaussianProposalInfo::get_covariance()
{
  return this->covariance;
}

double& GaussianProposalInfo::get_double_scale()
{
  return this->scale.get_constant();
}

Scale& GaussianProposalInfo::get_scale()
{
  return this->scale;
}

arma::mat& GaussianProposalInfo::get_chol()
{
  return this->chol;
}

arma::mat& GaussianProposalInfo::get_inv()
{
  return this->inv;
}

arma::mat& GaussianProposalInfo::get_inv_chol()
{
  return this->inv_chol;
}

double& GaussianProposalInfo::get_logdet()
{
  return this->logdet;
}
