#include <iterator>
#include "gaussian_independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

GaussianIndependentProposalKernel::GaussianIndependentProposalKernel()
  :IndependentProposalKernel()
{
}

GaussianIndependentProposalKernel::~GaussianIndependentProposalKernel()
{
}

// find mean and cov adaptively
GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in)
  :IndependentProposalKernel()
{
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

// find cov adaptively
GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<arma::colvec> &means_in)
  :IndependentProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(means_in[i]);
  }
}

// find mean adaptively
GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<arma::mat> &covariances_in)
  :IndependentProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
}

GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::string &variable_name_in,
                                                                     const double &mean_in,
                                                                     const double &sd_in)
:IndependentProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,sd_in);
}

GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::string &variable_name_in,
                                                                     const arma::colvec &mean_in,
                                                                     const arma::mat &covariance_in)
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,covariance_in);
}

GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<arma::colvec> &means_in,
                                                                     const std::vector<arma::mat> &covariances_in)
  :IndependentProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(means_in[i],covariances_in[i]);
  }
}

GaussianIndependentProposalKernel::GaussianIndependentProposalKernel(const GaussianIndependentProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void GaussianIndependentProposalKernel::operator=(const GaussianIndependentProposalKernel &another)
{
  if(this == &another)
    return;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* GaussianIndependentProposalKernel::duplicate() const
{
  return( new GaussianIndependentProposalKernel(*this));
}

ProposalKernel* GaussianIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new GaussianIndependentProposalKernel(*this));
}

IndependentProposalKernel* GaussianIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new GaussianIndependentProposalKernel(*this));
}

void GaussianIndependentProposalKernel::make_copy(const GaussianIndependentProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
}

void GaussianIndependentProposalKernel::set_mean(const std::string &variable,
                                                 const arma::colvec &mean_in)
{
  this->proposal_info[variable].set_mean(mean_in);
}

void GaussianIndependentProposalKernel::set_covariance(const std::string &variable,
                                                       const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat GaussianIndependentProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat GaussianIndependentProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}

arma::mat GaussianIndependentProposalKernel::get_covariance(const std::vector<std::string> &variables)
{
  size_t total_dimension = 0;
  std::vector<arma::mat> covariances;
  covariances.reserve(variables.size());
  for (auto i=variables.begin();
       i!=variables.end();
       ++i)
  {
    covariances.push_back(this->get_covariance(*i));
    total_dimension = total_dimension + covariances.back().n_rows;
  }
  
  arma::mat covariance(total_dimension,total_dimension);
  size_t current_dimension = 0;
  for (auto i=covariances.begin();
       i!=covariances.end();
       ++i)
  {
    covariance( arma::span(current_dimension,current_dimension+(*i).n_rows-1), arma::span(current_dimension,current_dimension+(*i).n_rows-1) ) = *i;
  }
  
  return covariance;
}

double GaussianIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_colvec(i->first),
                                            mean,
                                            (1.0/scale)*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
double GaussianIndependentProposalKernel::evaluate_independent_kernel(Variables* proposed_particle,
                                                                      const Parameters &conditioned_on_parameters) const
{
  return this->evaluate_independent_kernel(proposed_particle);
}
*/

double GaussianIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->evaluate_independent_kernel(proposed_particle);
}

Parameters GaussianIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    output[i->first] = rmvnorm(rng,
                               i->second.get_mean(),
                               sqrt(i->second.get_double_scale())*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters GaussianIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                   const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng);
}

Parameters GaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters GaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters GaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  output[variable] = rmvnorm(rng,
                             found->second.get_mean(),
                             sqrt(found->second.get_double_scale())*found->second.get_chol(),
                             true);
  return output;
}

Parameters GaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const std::string &variable,
                                                                             const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  output[variable] = rmvnorm(rng,
                             found->second.get_mean(),
                             sqrt(found->second.get_double_scale())*found->second.get_chol(),
                             true);
  return output;
}

arma::mat GaussianIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                         const Parameters &proposed_particle)
{
  Rcpp::stop("GaussianIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

arma::mat GaussianIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                   const Parameters &proposed_particle)
{
  Rcpp::stop("GaussianIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

void GaussianIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}
