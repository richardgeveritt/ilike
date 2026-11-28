#include <iterator>
#include "normal_source.h"
#include "nonlinear_gaussian_noise_function_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{
NonLinearGaussianNoiseFunctionProposalKernel::NonLinearGaussianNoiseFunctionProposalKernel()
:GaussianNoiseProposalKernel()
{
}

NonLinearGaussianNoiseFunctionProposalKernel::~NonLinearGaussianNoiseFunctionProposalKernel()
{
}

NonLinearGaussianNoiseFunctionProposalKernel::NonLinearGaussianNoiseFunctionProposalKernel(const std::string &variable_name_in,
                                                                                           std::shared_ptr<Transform> transform_in,
                                                                                           const GetMatrixPtr &covariance_in)
:GaussianNoiseProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalFunctionsInfo();
  this->proposal_info[variable_name_in].set_covariance(covariance_in);
  this->transform = transform_in;
}

NonLinearGaussianNoiseFunctionProposalKernel::NonLinearGaussianNoiseFunctionProposalKernel(const NonLinearGaussianNoiseFunctionProposalKernel &another)
:GaussianNoiseProposalKernel(another)
{
  this->make_copy(another);
}

void NonLinearGaussianNoiseFunctionProposalKernel::operator=(const NonLinearGaussianNoiseFunctionProposalKernel &another)
{
  if(this == &another)
    return;
  
  GaussianNoiseProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* NonLinearGaussianNoiseFunctionProposalKernel::duplicate() const
{
  return( new NonLinearGaussianNoiseFunctionProposalKernel(*this));
}

ProposalKernel* NonLinearGaussianNoiseFunctionProposalKernel::proposal_kernel_duplicate() const
{
  return( new NonLinearGaussianNoiseFunctionProposalKernel(*this));
}

GaussianNoiseProposalKernel* NonLinearGaussianNoiseFunctionProposalKernel::gaussian_noise_proposal_kernel_duplicate() const
{
  return( new NonLinearGaussianNoiseFunctionProposalKernel(*this));
}

void NonLinearGaussianNoiseFunctionProposalKernel::make_copy(const NonLinearGaussianNoiseFunctionProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  this->transform = another.transform;
}

double NonLinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                              const Particle &old_particle) const
{
  return this->specific_evaluate_kernel_from_parameters(proposed_particle.get_transformed_parameters(this),
                                                        old_particle.get_transformed_parameters(this));
}

double NonLinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
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
                                             i->second.get_inv(transformed),
                                             i->second.get_double_scale(),
                                             i->second.get_logdet(transformed));
  }
  return output;
}

double NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                        const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

// (dev/plan.md, PMMH P2) simulate, written into the target: the same random
// numbers and values, without a new Parameters for the output
void NonLinearGaussianNoiseFunctionProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                                 const Particle &particle,
                                                                 ProposalTarget &target) const
{
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    target.matrix(i->first) = source_rmvnorm_using_chol(rng,
                                                 transformed.get_colvec(i->first),
                                                 sqrt(scale)*i->second.get_chol(transformed));
  }
}

Parameters NonLinearGaussianNoiseFunctionProposalKernel::simulate(RandomNumberGenerator &rng,
                                                                  const Particle &particle) const
{
  Parameters output;
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    output[i->first] = source_rmvnorm_using_chol(rng,
                                          transformed.get_colvec(i->first),
                                          sqrt(scale)*i->second.get_chol(transformed));
  }
  return output;
}

Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                            const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                            const std::string &variable,
                                                                            const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  
  Parameters output;
  output[variable] = source_rmvnorm_using_chol(rng,
                                        transformed.get_colvec(variable),
                                        sqrt(found->second.get_double_scale())*found->second.get_chol(transformed));
  return output;
}

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                                 const Particle &proposed_particle,
                                                                                 const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                           const Particle &proposed_particle,
                                                                                           const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

void NonLinearGaussianNoiseFunctionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

std::vector<std::string> NonLinearGaussianNoiseFunctionProposalKernel::get_variables() const
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

GradientEstimatorOutput* NonLinearGaussianNoiseFunctionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> NonLinearGaussianNoiseFunctionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void NonLinearGaussianNoiseFunctionProposalKernel::set_index(Index* index_in)
{
}

void NonLinearGaussianNoiseFunctionProposalKernel::set_index_if_null(Index* index_in)
{
}

// (C269, PMMH P1)
void NonLinearGaussianNoiseFunctionProposalKernel::set_covariance_constant()
{
  for (auto i=this->proposal_info.begin(); i!=this->proposal_info.end(); ++i)
    i->second.set_covariance_constant();
}

// (C269)
void NonLinearGaussianNoiseFunctionProposalKernel::set_declared_reads(const std::vector<std::string> &reads)
{
  for (auto i=this->proposal_info.begin(); i!=this->proposal_info.end(); ++i)
    i->second.set_declared_reads(reads);
}

void NonLinearGaussianNoiseFunctionProposalKernel::prepare_for_parameters(const Parameters &parameters)
{
  for (auto i=this->proposal_info.begin(); i!=this->proposal_info.end(); ++i)
    i->second.prepare(parameters);
}

}
