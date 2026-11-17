#include "hmc_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{
HMCProposalKernel::HMCProposalKernel()
:ProposalKernel()
{
  
}

HMCProposalKernel::~HMCProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  std::vector<arma::mat> covariances_in;
  std::vector<std::string> covariance_names_in;
  covariances_in.reserve(variable_names_in.size());
  covariance_names_in.reserve(variable_names_in.size());
  for (size_t i=0; i<variable_names_in.size(); ++i)
  {
    covariances_in.push_back(arma::mat());
    covariance_names_in.push_back(variable_names_in[i] + "_cov");
  }
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                     const std::vector<arma::mat> &covariances_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  std::vector<std::string> covariance_names_in;
  covariance_names_in.reserve(variable_names_in.size());
  for (size_t i=0; i<variable_names_in.size(); ++i)
  {
    covariance_names_in.push_back(variable_names_in[i] + "_cov");
  }
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                     const std::vector<std::string> &covariance_names_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  std::vector<arma::mat> covariances_in;
  covariances_in.reserve(covariance_names_in.size());
  for (size_t i=0; i<covariance_names_in.size(); ++i)
  {
    covariances_in.push_back(arma::mat());
  }
}

HMCProposalKernel::HMCProposalKernel(const std::string &variable_name_in,
                                     const arma::mat &covariance_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const std::string &variable_name_in,
                                     const double &sd_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                     const std::vector<std::string> &covariance_names_in,
                                     const std::vector<arma::mat> &covariances_in,
                                     GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const HMCProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void HMCProposalKernel::operator=(const HMCProposalKernel &another)
{
  if(this == &another)
    return;
  
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* HMCProposalKernel::duplicate() const
{
  return( new HMCProposalKernel(*this));
}

ProposalKernel* HMCProposalKernel::proposal_kernel_duplicate() const
{
  return( new HMCProposalKernel(*this));
}

void HMCProposalKernel::make_copy(const HMCProposalKernel &another)
{
  this->variable_names = another.variable_names;
  this->covariance_names = another.covariance_names;
}

double HMCProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                   const Particle &old_particle) const
{
  double output = 0.0;
  return output;
}

double HMCProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                             const Particle &old_particle) const
{
  // needs changing
  double output = 0.0;
  return output;
}

Parameters HMCProposalKernel::simulate(RandomNumberGenerator &rng,
                                       const Particle &particle) const
{
  Parameters output;
  return output;
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const Particle &particle) const
{
  // needs changing
  Parameters output;
  return output;
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const std::string &variable,
                                                 const Particle &particle) const
{
  // needs changing
  Parameters output;
  return output;
}

arma::mat HMCProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                      const Particle &proposed_particle,
                                                      const Particle &old_particle)
{
  Rcpp::stop("HMCProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat HMCProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                const Particle &proposed_particle,
                                                                const Particle &old_particle)
{
  Rcpp::stop("HMCProposalKernel::specific_gradient_of_log - not written yet.");
}

void HMCProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

GradientEstimatorOutput* HMCProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> HMCProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void HMCProposalKernel::set_index(Index* index_in)
{
}

void HMCProposalKernel::set_index_if_null(Index* index_in)
{
}

bool HMCProposalKernel::can_be_evaluated() const
{
  return true;
}

void HMCProposalKernel::set_data(Data* data_in)
{
  
}
}
