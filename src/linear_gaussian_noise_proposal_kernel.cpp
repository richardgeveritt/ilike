#include <iterator>
#include "normal_source.h"
#include "linear_gaussian_noise_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{
LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel()
:GaussianNoiseProposalKernel()
{
}

LinearGaussianNoiseProposalKernel::~LinearGaussianNoiseProposalKernel()
{
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                     const std::vector<std::string> &conditioned_on_variable_names_in,
                                                                     const std::vector<arma::mat> &As_in)
:GaussianNoiseProposalKernel()
{
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
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
  this->conditioned_on_variable_names = conditioned_on_variable_names_in;
}

LinearGaussianNoiseProposalKernel::LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                                                     const std::string &conditioned_on_variable_name_in,
                                                                     const arma::mat &A_in,
                                                                     const arma::mat &covariance_in)
:GaussianNoiseProposalKernel()
{
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
}

double LinearGaussianNoiseProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                   const Particle &old_particle) const
{
  return this->specific_evaluate_kernel_from_parameters(proposed_particle.get_transformed_parameters(this),
                                                        old_particle.get_transformed_parameters(this));
}

double LinearGaussianNoiseProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                                   const Parameters &old_parameters) const
{
  double output = 0.0;

  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_A()*old_parameters.get_colvec(conditioned_on_variable_names);
    output = output + dmvnorm_scaled_precomp(proposed_parameters.get_colvec(i->first),
                                             mean,
                                             i->second.get_inv(),
                                             i->second.get_double_scale(),
                                             i->second.get_logdet());
  }
  return output;
}

double LinearGaussianNoiseProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                             const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

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

// simulate, written into the target (same random numbers and values)
void LinearGaussianNoiseProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                      const Particle &particle,
                                                      ProposalTarget &target) const
{
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    target.matrix(i->first) = source_rmvnorm_using_lower_chol(rng,
                                                       i->second.get_A()*particle.get_transformed_parameters(this).get_colvec(i->first),
                                                       i->second.get_chol_lower(), sqrt(scale));
  }
}

Parameters LinearGaussianNoiseProposalKernel::simulate(RandomNumberGenerator &rng,
                                                       const Particle &particle) const
{
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    double scale = i->second.get_double_scale();
    output[i->first] = source_rmvnorm_using_lower_chol(rng,
                                          i->second.get_A()*particle.get_transformed_parameters(this).get_colvec(i->first),
                                          i->second.get_chol_lower(), sqrt(scale));
  }
  return output;
}

Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                 const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

Parameters LinearGaussianNoiseProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                 const std::string &variable,
                                                                 const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  output[variable] = source_rmvnorm_using_lower_chol(rng,
                                        found->second.get_A()*particle.get_transformed_parameters(this).get_colvec(variable),
                                        found->second.get_chol_lower(), sqrt(found->second.get_double_scale()));
  return output;
}

arma::mat LinearGaussianNoiseProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                      const Particle &proposed_particle,
                                                                      const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat LinearGaussianNoiseProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                const Particle &proposed_particle,
                                                                                const Particle &old_particle)
{
  Rcpp::stop("LinearGaussianNoiseProposalKernel::specific_gradient_of_log - not written yet.");
}

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

// (dev/plan.md 3.21) a Gaussian density is largest at its mean, whatever the
// mean (here a function of the state): the density at the mean, for each
// variable, with the covariances fixed
double LinearGaussianNoiseProposalKernel::log_density_maximum(const Parameters &parameters) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin(); i!=this->proposal_info.end(); ++i)
  {
    arma::colvec at_mean(i->second.get_inv().n_rows, arma::fill::zeros);
    output = output + dmvnorm_scaled_precomp(at_mean, at_mean, i->second.get_inv(),
                                             i->second.get_double_scale(), i->second.get_logdet());
  }
  return output;
}
}
