#include "gaussian_proposal_functions_info.h"
#include "utils.h"

GaussianProposalFunctionsInfo::GaussianProposalFunctionsInfo()
: scale(1.0)
{
}

GaussianProposalFunctionsInfo::~GaussianProposalFunctionsInfo()
{
  
}

//Copy constructor for the GaussianProposalFunctionsInfo class.
GaussianProposalFunctionsInfo::GaussianProposalFunctionsInfo(const GaussianProposalFunctionsInfo &another)
{
  this->make_copy(another);
}

GaussianProposalFunctionsInfo& GaussianProposalFunctionsInfo::operator=(const GaussianProposalFunctionsInfo &another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(another);
  return *this;
}

GaussianProposalFunctionsInfo* GaussianProposalFunctionsInfo::duplicate() const
{
  return( new GaussianProposalFunctionsInfo(*this));
}

void GaussianProposalFunctionsInfo::make_copy(const GaussianProposalFunctionsInfo &another)
{
  this->covariance = another.covariance;
  this->A = another.A;
  this->mean = another.mean;
  this->scale = another.scale;
  this->chol = another.chol;
}

//Move constructor for the GaussianProposalFunctionsInfo class.
GaussianProposalFunctionsInfo::GaussianProposalFunctionsInfo(GaussianProposalFunctionsInfo &&another)
: mean(NULL), covariance(NULL), scale(NULL)
{
  this->make_copy(std::move(another));
}

GaussianProposalFunctionsInfo& GaussianProposalFunctionsInfo::operator=(GaussianProposalFunctionsInfo &&another)
{
  if(this == &another){ //if a==a
    return *this;
  }
  
  this->make_copy(std::move(another));
  return *this;
}

void GaussianProposalFunctionsInfo::make_copy(GaussianProposalFunctionsInfo &&another)
{
  this->covariance = another.covariance;
  this->A = another.A;
  this->mean = another.mean;
  this->scale = std::move(another.scale);
  
  this->chol = std::move(another.chol);
  
  another.covariance = NULL;
  another.A = NULL;
  another.mean = NULL;
  another.scale = Scale();
  
  another.chol = arma::mat();
}

void GaussianProposalFunctionsInfo::set_mean(const GetVectorPtr &mean_in)
{
  this->mean = mean_in;
}

void GaussianProposalFunctionsInfo::set_covariance(const GetMatrixPtr &covariance_in)
{
  this->covariance = covariance_in;
}

void GaussianProposalFunctionsInfo::set_A(const GetMatrixPtr &A_in)
{
  this->A = A_in;
}

arma::colvec GaussianProposalFunctionsInfo::get_mean(const Parameters &parameters) const
{
  return this->mean(parameters);
}

arma::mat GaussianProposalFunctionsInfo::get_covariance(const Parameters &parameters) const
{
  return this->covariance(parameters);
}

arma::mat GaussianProposalFunctionsInfo::get_A(const Parameters &parameters) const
{
  return this->A(parameters);
}

double GaussianProposalFunctionsInfo::get_double_scale() const
{
  return this->scale.get_constant();
}

Scale GaussianProposalFunctionsInfo::get_scale() const
{
  return this->scale;
}

double& GaussianProposalFunctionsInfo::get_double_scale()
{
  return this->scale.get_constant();
}

Scale& GaussianProposalFunctionsInfo::get_scale()
{
  return this->scale;
}

arma::mat GaussianProposalFunctionsInfo::get_chol(const Parameters &parameters) const
{
  return arma::chol(this->covariance(parameters));
}

arma::mat GaussianProposalFunctionsInfo::get_inv(const Parameters &parameters) const
{
  return arma::inv_sympd(this->covariance(parameters));
}

double GaussianProposalFunctionsInfo::get_logdet(const Parameters &parameters) const
{
  return arma::log_det_sympd(this->covariance(parameters));
}
