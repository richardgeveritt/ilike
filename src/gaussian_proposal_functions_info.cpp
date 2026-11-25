#include "gaussian_proposal_functions_info.h"
#include "utils.h"

namespace ilike
{
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
  this->covariance_constant = another.covariance_constant;
  this->constant_chol = another.constant_chol;
  this->constant_inv = another.constant_inv;
  this->constant_logdet = another.constant_logdet;
  this->declared_reads = another.declared_reads;
  this->cache_valid = another.cache_valid;
  this->cached_values = another.cached_values;
  this->cached_A = another.cached_A;
  this->cached_A_set = another.cached_A_set;
  this->cached_chol = another.cached_chol;
  this->cached_inv = another.cached_inv;
  this->cached_logdet = another.cached_logdet;
  this->cached_covariance_set = another.cached_covariance_set;
}

//Move constructor for the GaussianProposalFunctionsInfo class.
GaussianProposalFunctionsInfo::GaussianProposalFunctionsInfo(GaussianProposalFunctionsInfo &&another)
: mean(NULL), covariance(NULL), scale()
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
  this->covariance_constant = another.covariance_constant;
  this->constant_chol = std::move(another.constant_chol);
  this->constant_inv = std::move(another.constant_inv);
  this->constant_logdet = another.constant_logdet;
  this->declared_reads = std::move(another.declared_reads);
  this->cache_valid = another.cache_valid;
  this->cached_values = std::move(another.cached_values);
  this->cached_A = std::move(another.cached_A);
  this->cached_A_set = another.cached_A_set;
  this->cached_chol = std::move(another.cached_chol);
  this->cached_inv = std::move(another.cached_inv);
  this->cached_logdet = another.cached_logdet;
  this->cached_covariance_set = another.cached_covariance_set;
  
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

void GaussianProposalFunctionsInfo::set_covariance_constant()
{
  // (the covariance function returns the constant whatever its parameters)
  arma::mat value = this->covariance(Parameters());
  this->constant_chol = arma::chol(value);
  this->constant_inv = arma::inv_sympd(value);
  this->constant_logdet = arma::log_det_sympd(value);
  this->covariance_constant = true;
}

void GaussianProposalFunctionsInfo::set_declared_reads(const std::vector<std::string> &reads_in)
{
  this->declared_reads = reads_in;
  this->cache_valid = false;
}

void GaussianProposalFunctionsInfo::prepare(const Parameters &parameters)
{
  if (this->declared_reads.empty())
    return;
  this->cache_valid = false;
  // only the declared variables: a function reading another one stops here
  Parameters declared;
  std::vector<arma::mat> values;
  values.reserve(this->declared_reads.size());
  for (const std::string &name : this->declared_reads)
  {
    const arma::mat* value = parameters.find_vector(name);
    if (value == NULL)
      return; // (the declared variables are not known here: no cache)
    declared[name] = *value;
    values.push_back(*value);
  }
  try
  {
    this->cached_A_set = false;
    if (this->A != NULL)
    {
      this->cached_A = this->A(declared);
      this->cached_A_set = true;
    }
    this->cached_covariance_set = false;
    if (!this->covariance_constant && (this->covariance != NULL))
    {
      arma::mat value = this->covariance(declared);
      this->cached_chol = arma::chol(value);
      this->cached_inv = arma::inv_sympd(value);
      this->cached_logdet = arma::log_det_sympd(value);
      this->cached_covariance_set = true;
    }
  }
  catch (std::exception &e)
  {
    Rcpp::stop(std::string("a transition's reads = [...] declares what its matrix and covariance functions read, but "
                           "evaluating them on only those variables failed (do they read another variable?): ") + e.what());
  }
  this->cached_values = std::move(values);
  this->cache_valid = true;
}

bool GaussianProposalFunctionsInfo::cache_matches(const Parameters &parameters) const
{
  if (!this->cache_valid)
    return false;
  for (size_t k = 0; k < this->declared_reads.size(); ++k)
  {
    const arma::mat* value = parameters.find_vector(this->declared_reads[k]);
    if (value == NULL)
      return false;
    const arma::mat &kept = this->cached_values[k];
    if ((value->n_rows != kept.n_rows) || (value->n_cols != kept.n_cols) || arma::any(arma::vectorise(*value != kept)))
      return false;
  }
  return true;
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
  if (this->cached_A_set && this->cache_matches(parameters))
    return this->cached_A;
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
  if (this->covariance_constant)
    return this->constant_chol;
  if (this->cached_covariance_set && this->cache_matches(parameters))
    return this->cached_chol;
  return arma::chol(this->covariance(parameters));
}

arma::mat GaussianProposalFunctionsInfo::get_inv(const Parameters &parameters) const
{
  if (this->covariance_constant)
    return this->constant_inv;
  if (this->cached_covariance_set && this->cache_matches(parameters))
    return this->cached_inv;
  return arma::inv_sympd(this->covariance(parameters));
}

double GaussianProposalFunctionsInfo::get_logdet(const Parameters &parameters) const
{
  if (this->covariance_constant)
    return this->constant_logdet;
  if (this->cached_covariance_set && this->cache_matches(parameters))
    return this->cached_logdet;
  return arma::log_det_sympd(this->covariance(parameters));
}

}
