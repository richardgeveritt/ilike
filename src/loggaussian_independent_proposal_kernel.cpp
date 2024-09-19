#include <iterator>
#include "loggaussian_independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

namespace ilike
{
LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel()
:IndependentProposalKernel()
{
}

LogGaussianIndependentProposalKernel::~LogGaussianIndependentProposalKernel()
{
}

// find mean and cov adaptively
LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in)
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
LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
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
LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
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

LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::string &variable_name_in,
                                                                           const double &mean_in,
                                                                           const double &sd_in)
:IndependentProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,sd_in);
}

LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::string &variable_name_in,
                                                                           const arma::colvec &mean_in,
                                                                           const arma::mat &sd_in)
:IndependentProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(mean_in,sd_in);
}

LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
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

LogGaussianIndependentProposalKernel::LogGaussianIndependentProposalKernel(const LogGaussianIndependentProposalKernel &another)
:IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void LogGaussianIndependentProposalKernel::operator=(const LogGaussianIndependentProposalKernel &another)
{
  if(this == &another)
    return;
  
  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* LogGaussianIndependentProposalKernel::duplicate() const
{
  return( new LogGaussianIndependentProposalKernel(*this));
}

ProposalKernel* LogGaussianIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new LogGaussianIndependentProposalKernel(*this));
}

IndependentProposalKernel* LogGaussianIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new LogGaussianIndependentProposalKernel(*this));
}

void LogGaussianIndependentProposalKernel::make_copy(const LogGaussianIndependentProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
}

void LogGaussianIndependentProposalKernel::set_mean(const std::string &variable,
                                                    const arma::colvec &mean_in)
{
  this->proposal_info[variable].set_mean(mean_in);
}

void LogGaussianIndependentProposalKernel::set_covariance(const std::string &variable,
                                                          const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat LogGaussianIndependentProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat LogGaussianIndependentProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}

arma::mat LogGaussianIndependentProposalKernel::get_covariance(const std::vector<std::string> &variables)
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

double LogGaussianIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvlnorm_using_precomp(proposed_particle.get_colvec(i->first),
                                             mean,
                                             (1.0/scale)*i->second.get_inv(),
                                             dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
 double LogGaussianIndependentProposalKernel::evaluate_independent_kernel(Variables* proposed_particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->evaluate_independent_kernel(proposed_particle);
 }
 */

double LogGaussianIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->evaluate_independent_kernel(proposed_particle);
}

Parameters LogGaussianIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    output[i->first] = rmvlnorm_using_chol(rng,
                                           i->second.get_mean(),
                                           sqrt(i->second.get_double_scale())*i->second.get_chol());
  }
  return output;
}

Parameters LogGaussianIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                      const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng);
}

Parameters LogGaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters LogGaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters LogGaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  output[variable] = rmvlnorm_using_chol(rng,
                                         found->second.get_mean(),
                                         sqrt(found->second.get_double_scale())*found->second.get_chol());
  return output;
}

Parameters LogGaussianIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const std::string &variable,
                                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  output[variable] = rmvlnorm_using_chol(rng,
                                         found->second.get_mean(),
                                         sqrt(found->second.get_double_scale())*found->second.get_chol());
  return output;
}

arma::mat LogGaussianIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                            const Parameters &proposed_particle)
{
  Rcpp::stop("LogGaussianIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

arma::mat LogGaussianIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                      const Parameters &proposed_particle)
{
  Rcpp::stop("LogGaussianIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

void LogGaussianIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

GradientEstimatorOutput* LogGaussianIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> LogGaussianIndependentProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void LogGaussianIndependentProposalKernel::set_index(Index* index_in)
{
}

void LogGaussianIndependentProposalKernel::set_index_if_null(Index* index_in)
{
}

bool LogGaussianIndependentProposalKernel::can_be_evaluated() const
{
  return true;
}

void LogGaussianIndependentProposalKernel::set_data(Data* data_in)
{
  
}
}
