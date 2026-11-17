#include <iterator>
#include "linear_gaussian_noise_function_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{
LinearGaussianNoiseFunctionProposalKernel::LinearGaussianNoiseFunctionProposalKernel()
:GaussianNoiseProposalKernel()
{
}

LinearGaussianNoiseFunctionProposalKernel::~LinearGaussianNoiseFunctionProposalKernel()
{
}

LinearGaussianNoiseFunctionProposalKernel::LinearGaussianNoiseFunctionProposalKernel(const std::string &variable_name_in,
                                                                                     const std::string &conditioned_on_variable_name_in,
                                                                                     const GetMatrixPtr &A_in,
                                                                                     const GetMatrixPtr &covariance_in)
:GaussianNoiseProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalFunctionsInfo();
  this->proposal_info[variable_name_in].set_covariance(covariance_in);
  this->proposal_info[variable_name_in].set_A(A_in);
  this->conditioned_on_variable_names.push_back(conditioned_on_variable_name_in);
}

LinearGaussianNoiseFunctionProposalKernel::LinearGaussianNoiseFunctionProposalKernel(const std::string &variable_name_in,
                                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                                     const GetMatrixPtr &A_in,
                                                                                     const GetMatrixPtr &covariance_in)
:GaussianNoiseProposalKernel()
{
  this->proposal_info[variable_name_in] = GaussianProposalFunctionsInfo();
  this->proposal_info[variable_name_in].set_covariance(covariance_in);
  this->proposal_info[variable_name_in].set_A(A_in);
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseFunctionProposalKernel::LinearGaussianNoiseFunctionProposalKernel(const LinearGaussianNoiseFunctionProposalKernel &another)
:GaussianNoiseProposalKernel(another)
{
  this->make_copy(another);
}

void LinearGaussianNoiseFunctionProposalKernel::operator=(const LinearGaussianNoiseFunctionProposalKernel &another)
{
  if(this == &another)
    return;
  
  GaussianNoiseProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* LinearGaussianNoiseFunctionProposalKernel::duplicate() const
{
  return( new LinearGaussianNoiseFunctionProposalKernel(*this));
}

ProposalKernel* LinearGaussianNoiseFunctionProposalKernel::proposal_kernel_duplicate() const
{
  return( new LinearGaussianNoiseFunctionProposalKernel(*this));
}

GaussianNoiseProposalKernel* LinearGaussianNoiseFunctionProposalKernel::gaussian_noise_proposal_kernel_duplicate() const
{
  return( new LinearGaussianNoiseFunctionProposalKernel(*this));
}

void LinearGaussianNoiseFunctionProposalKernel::make_copy(const LinearGaussianNoiseFunctionProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  this->conditioned_on_variable_names = another.conditioned_on_variable_names;
}

double LinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                           const Particle &old_particle) const
{
  return this->specific_evaluate_kernel_from_parameters(proposed_particle.get_transformed_parameters(this),
                                                        old_particle.get_transformed_parameters(this));
}

double LinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                                           const Parameters &old_parameters) const
{
  double output = 0.0;

  const Parameters &current_parameters = old_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_A(current_parameters)*current_parameters.get_colvec(this->conditioned_on_variable_names);
    output = output + dmvnorm_scaled_precomp(proposed_parameters.get_colvec(i->first),
                                             mean,
                                             i->second.get_inv(current_parameters),
                                             i->second.get_double_scale(),
                                             i->second.get_logdet(current_parameters));
  }
  return output;
}

double LinearGaussianNoiseFunctionProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                     const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters LinearGaussianNoiseFunctionProposalKernel::simulate(RandomNumberGenerator &rng,
                                                               const Particle &particle) const
{
  Parameters output;
  Parameters current_parameters = particle.get_transformed_parameters(this);
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    output[i->first] = rmvnorm_using_chol(rng,
                                          i->second.get_A(current_parameters)*current_parameters.get_colvec(i->first),
                                          sqrt(scale)*i->second.get_chol(current_parameters));
  }
  return output;
}

Parameters LinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                         const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

Parameters LinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                         const std::string &variable,
                                                                         const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  Parameters current_parameters = particle.get_transformed_parameters(this);
  output[variable] = rmvnorm_using_chol(rng,
                                        found->second.get_A(current_parameters)*current_parameters.get_colvec(variable),
                                        sqrt(found->second.get_double_scale())*found->second.get_chol(current_parameters));
  return output;
}

arma::mat LinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                              const Particle &proposed_particle,
                                                                              const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat LinearGaussianNoiseFunctionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                        const Particle &proposed_particle,
                                                                                        const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

void LinearGaussianNoiseFunctionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

std::vector<std::string> LinearGaussianNoiseFunctionProposalKernel::get_variables() const
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

GradientEstimatorOutput* LinearGaussianNoiseFunctionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> LinearGaussianNoiseFunctionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void LinearGaussianNoiseFunctionProposalKernel::set_index(Index* index_in)
{
}

void LinearGaussianNoiseFunctionProposalKernel::set_index_if_null(Index* index_in)
{
}
}
