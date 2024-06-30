#include <iterator>
#include "linear_gaussian_noise_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel()
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
}

LinearGaussianNoiseProposalKernel::~LinearGaussianNoiseProposalKernel()
{
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                     const std::vector<arma::mat> &As_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo();
    this->proposal_info[variable_names_in[i]].set_A(As_in[i]);
  }
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                     const std::vector<arma::mat> &As_in,
                                                                     const std::vector<arma::mat> &covariances_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
    this->proposal_info[variable_names_in[i]].set_A(As_in[i]);
  }
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                     const double &sd_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                     const std::string &conditioned_on_variable_name_in,
                                                                     const arma::mat &A_in,
                                                                     const arma::mat &covariance_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  this->proposal_info[variable_name_in].set_A(A_in);
  this->conditioned_on_variable_names.push_back(conditioned_on_variable_name_in);
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                     const arma::mat &A_in,
                                                                     const arma::mat &covariance_in)
:GaussianNoiseProposalKernel()
{
  //this->unused_variables_kept = true;
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
  this->proposal_info[variable_name_in].set_A(A_in);
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const LinearGaussianNoiseProposalKernel &another)
:GaussianNoiseProposalKernel(another)
{
  this->make_copy(another);
}

void LinearGaussianNoiseProposalKernel::operator=(const LinearGaussianNoiseProposalKernel &another)
{
  if(this == &another)
    return;
  
  GaussianNoiseProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* LinearGaussianNoiseProposalKernel::duplicate() const
{
  return( new LinearGaussianNoiseProposalKernel(*this));
}

ProposalKernel* LinearGaussianNoiseProposalKernel::proposal_kernel_duplicate() const
{
  return( new LinearGaussianNoiseProposalKernel(*this));
}

GaussianNoiseProposalKernel* LinearGaussianNoiseProposalKernel::gaussian_noise_proposal_kernel_duplicate() const
{
  return( new LinearGaussianNoiseProposalKernel(*this));
}

void LinearGaussianNoiseProposalKernel::make_copy(const LinearGaussianNoiseProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  this->conditioned_on_variable_names = another.conditioned_on_variable_names;
  //this->unused_variables_kept = another.unused_variables_kept;
}

double LinearGaussianNoiseProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                   const Particle &old_particle) const
{
  double output = 0.0;

  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_A()*old_particle.get_transformed_parameters(this).get_colvec(conditioned_on_variable_names);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(i->first),
                                            mean,
                                            (1.0/sqrt(scale))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
 double LinearGaussianNoiseProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->specific_evaluate_kernel(proposed_particle, old_particle);
 }
 */

double LinearGaussianNoiseProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                             const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
 double LinearGaussianNoiseProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->specific_evaluate_kernel(proposed_particle, old_particle);
 }
 */

void LinearGaussianNoiseProposalKernel::set_covariance(const std::string &variable,
                                                       const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat LinearGaussianNoiseProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat LinearGaussianNoiseProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}

Parameters LinearGaussianNoiseProposalKernel::simulate(RandomNumberGenerator &rng,
                                                       const Particle &particle) const
{
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    //arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm_using_chol(rng,
                                          i->second.get_A()*particle.get_transformed_parameters(this).get_colvec(i->first),
                                          sqrt(scale)*i->second.get_chol());
  }
  return output;
}

/*
 Parameters LinearGaussianNoiseProposalKernel::simulate(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->simulate(rng, particle);
 }
 */

Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
 Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->simulate(rng, particle);
 }
 */

Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const std::string &variable,
                                                                const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  output[variable] = rmvnorm_using_chol(rng,
                                        found->second.get_A()*particle.get_transformed_parameters(this).get_colvec(variable),
                                        sqrt(found->second.get_double_scale())*found->second.get_chol());
  return output;
}

/*
 Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
 const std::string &variable,
 Particle &particle,
 const Parameters &conditioned_on_parameters) const
 {
 // no difference since size of data set does not impact on proposal
 return this->subsample_simulate(rng, variable, particle);
 }
 */

arma::mat LinearGaussianNoiseProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                     const Particle &proposed_particle,
                                                                     const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
 arma::mat LinearGaussianNoiseProposalKernel::specific_gradient_of_log(const std::string &variable,
 Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
 }
 */

arma::mat LinearGaussianNoiseProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                               const Particle &proposed_particle,
                                                                               const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
 arma::mat LinearGaussianNoiseProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
 Particle &proposed_particle,
 Particle &old_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
 }
 */

void LinearGaussianNoiseProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

std::vector<std::string> LinearGaussianNoiseProposalKernel::get_variables() const
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

GradientEstimatorOutput* LinearGaussianNoiseProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> LinearGaussianNoiseProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void LinearGaussianNoiseProposalKernel::set_index(Index* index_in)
{
}

void LinearGaussianNoiseProposalKernel::set_index_if_null(Index* index_in)
{
}
