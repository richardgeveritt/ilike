#include <iterator>
#include "nonlinear_gaussian_noise_function_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

NonLinearGaussianNoiseFunctionProposalKernel::NonLinearGaussianNoiseFunctionProposalKernel()
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
}

NonLinearGaussianNoiseFunctionProposalKernel::~NonLinearGaussianNoiseFunctionProposalKernel()
{
}

NonLinearGaussianNoiseFunctionProposalKernel::NonLinearGaussianNoiseFunctionProposalKernel(const std::string &variable_name_in,
                                                                                           std::shared_ptr<Transform> transform_in,
                                                                                           const const GetMatrixPtr &covariance_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  //this->unused_variables_kept = true;
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
  //this->unused_variables_kept = another.unused_variables_kept;
}

double NonLinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                              const Particle &old_particle) const
{
  double output = 0.0;
  Parameters transformed = this->transform->transform(old_particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = transformed.get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(i->first),
                                            mean,
                                            (1.0/sqrt(scale))*i->second.get_inv(transformed),
                                            dim*log(scale)+i->second.get_logdet(transformed));
  }
  return output;
}

/*
 double NonLinearGaussianNoiseFunctionProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->specific_evaluate_kernel(proposed_particle, old_particle);
 }
 */

double NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
 double NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->specific_evaluate_kernel(proposed_particle, old_particle);
 }
 */

/*
void NonLinearGaussianNoiseFunctionProposalKernel::set_covariance(const std::string &variable,
                                                          const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}
*/

Parameters NonLinearGaussianNoiseFunctionProposalKernel::simulate(RandomNumberGenerator &rng,
                                                                  const Particle &particle) const
{
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    //arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm_using_chol(rng,
                                          transformed.get_colvec(i->first),
                                          sqrt(scale)*i->second.get_chol(transformed));
  }
  return output;
}

/*
 Parameters NonLinearGaussianNoiseFunctionProposalKernel::simulate(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->simulate(rng, particle);
 }
 */

Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
 Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->simulate(rng, particle);
 }
 */

Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                    const std::string &variable,
                                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters transformed = this->transform->transform(particle.get_transformed_parameters(this));
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  output[variable] = rmvnorm_using_chol(rng,
                                        transformed.get_colvec(variable),
                                        sqrt(found->second.get_double_scale())*found->second.get_chol(transformed));
  return output;
}

/*
 Parameters NonLinearGaussianNoiseFunctionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
 const std::string &variable,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->subsample_simulate(rng, variable, particle);
 }
 */

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                         const Particle &proposed_particle,
                                                                         const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
 arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log(const std::string &variable,
 Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
 }
 */

arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                   const Particle &proposed_particle,
                                                                                   const Particle &old_particle)
{
  Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
 arma::mat NonLinearGaussianNoiseFunctionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
 Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("NonLinearGaussianNoiseFunctionProposalKernel::specific_gradient_of_log - not written yet.");
 }
 */

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
