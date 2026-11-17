#include <iterator>
#include "nonlinear_gaussian_noise_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{
NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel()
:GaussianNoiseProposalKernel()
{
}

NonLinearGaussianNoiseProposalKernel::~NonLinearGaussianNoiseProposalKernel()
{
}

NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                           std::shared_ptr<Transform> transform_in)
:GaussianNoiseProposalKernel()
{
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
  this->transform = transform_in;
}

NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                           std::shared_ptr<Transform> transform_in,
                                                                           const std::vector<arma::mat> &covariances_in)
:GaussianNoiseProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
  this->transform = transform_in;
}

NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                           std::shared_ptr<Transform> transform_in,
                                                                           const double &sd_in)
:GaussianNoiseProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
  this->transform = transform_in;
}

NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                           std::shared_ptr<Transform> transform_in,
                                                                           const arma::mat &covariance_in)
:GaussianNoiseProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  this->transform = transform_in;
}

NonLinearGaussianNoiseProposalKernel::NonLinearGaussianNoiseProposalKernel(const NonLinearGaussianNoiseProposalKernel &another)
:GaussianNoiseProposalKernel(another)
{
  this->make_copy(another);
}

void NonLinearGaussianNoiseProposalKernel::operator=(const NonLinearGaussianNoiseProposalKernel &another)
{
  if(this == &another)
    return;
  
  GaussianNoiseProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* NonLinearGaussianNoiseProposalKernel::duplicate() const
{
  return( new NonLinearGaussianNoiseProposalKernel(*this));
}

ProposalKernel* NonLinearGaussianNoiseProposalKernel::proposal_kernel_duplicate() const
{
  return( new NonLinearGaussianNoiseProposalKernel(*this));
}

GaussianNoiseProposalKernel* NonLinearGaussianNoiseProposalKernel::gaussian_noise_proposal_kernel_duplicate() const
{
  return( new NonLinearGaussianNoiseProposalKernel(*this));
}

void NonLinearGaussianNoiseProposalKernel::make_copy(const NonLinearGaussianNoiseProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  this->transform = another.transform;
}

double NonLinearGaussianNoiseProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                      const Particle &old_particle) const
{
  return this->specific_evaluate_kernel_from_parameters(proposed_particle.get_transformed_parameters(this),
                                                        old_particle.get_transformed_parameters(this));
}

double NonLinearGaussianNoiseProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                                      const Parameters &old_parameters) const
{
  double output = 0.0;
  Parameters transformed = this->transform->transform(old_parameters);
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    output = output + dmvnorm_scaled_precomp(proposed_parameters.get_colvec(i->first),
                                             transformed.get_colvec(i->first),
                                             i->second.get_inv(),
                                             i->second.get_double_scale(),
                                             i->second.get_logdet());
  }
  return output;
}

double NonLinearGaussianNoiseProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

void NonLinearGaussianNoiseProposalKernel::set_covariance(const std::string &variable,
                                                          const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat NonLinearGaussianNoiseProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat NonLinearGaussianNoiseProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}

// simulate, written into the target (same random numbers and values)
void NonLinearGaussianNoiseProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                         const Particle &particle,
                                                         ProposalTarget &target) const
{
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    target.matrix(i->first) = rmvnorm_using_lower_chol(rng,
                                                       transformed.get_colvec(i->first),
                                                       i->second.get_chol_lower(), sqrt(scale));
  }
}

Parameters NonLinearGaussianNoiseProposalKernel::simulate(RandomNumberGenerator &rng,
                                                          const Particle &particle) const
{
  Parameters output;
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    output[i->first] = rmvnorm_using_lower_chol(rng,
                                          transformed.get_colvec(i->first),
                                          i->second.get_chol_lower(), sqrt(scale));
  }
  return output;
}

Parameters NonLinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

Parameters NonLinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                    const std::string &variable,
                                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  
  Parameters output;
  output[variable] = rmvnorm_using_lower_chol(rng,
                                        transformed.get_colvec(variable),
                                        found->second.get_chol_lower(), sqrt(found->second.get_double_scale()));
  return output;
}

arma::mat NonLinearGaussianNoiseProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                         const Particle &proposed_particle,
                                                                         const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat NonLinearGaussianNoiseProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                   const Particle &proposed_particle,
                                                                                   const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

void NonLinearGaussianNoiseProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

std::vector<std::string> NonLinearGaussianNoiseProposalKernel::get_variables() const
{
  std::vector<std::string> variables;
  variables.reserve(this->proposal_info.size());
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    variables.push_back(i->first);
  }
  return variables;
}

GradientEstimatorOutput* NonLinearGaussianNoiseProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> NonLinearGaussianNoiseProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void NonLinearGaussianNoiseProposalKernel::set_index(Index* index_in)
{
}

void NonLinearGaussianNoiseProposalKernel::set_index_if_null(Index* index_in)
{
}
}
