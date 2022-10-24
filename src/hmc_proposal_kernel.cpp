#include "hmc_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

HMCProposalKernel::HMCProposalKernel()
  :ProposalKernel()
{
  
}

HMCProposalKernel::~HMCProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in)
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
                                     const std::vector<arma::mat> &covariances_in)
{
  std::vector<std::string> covariance_names_in;
  covariance_names_in.reserve(variable_names_in.size());
  for (size_t i=0; i<variable_names_in.size(); ++i)
  {
    covariance_names_in.push_back(variable_names_in[i] + "_cov");
  }
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                   const std::vector<std::string> &covariance_names_in)
{
  std::vector<arma::mat> covariances_in;
  covariances_in.reserve(covariance_names_in.size());
  for (size_t i=0; i<covariance_names_in.size(); ++i)
  {
    covariances_in.push_back(arma::mat());
  }
}

HMCProposalKernel::HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                   const std::vector<std::string> &covariance_names_in,
                                                                   const std::vector<arma::mat> &covariances_in)
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

double HMCProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                          Particle &old_particle) const
{
  double output = 0.0;
  /*
  for (size_t i=0; i<this->covariance_names.size(); ++i)
  {
    output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
                              old_particle.parameters.get_vector(this->variable_names[i]),
                              this->proposal_parameters[this->covariance_names[i]]);
  }
  */
  return output;
}

double HMCProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                          Particle &old_particle,
                                          const Parameters &conditioned_on_parameters) const
{
  double output = 0.0;
  /*
   for (size_t i=0; i<this->covariance_names.size(); ++i)
   {
   output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
   old_particle.parameters.get_vector(this->variable_names[i]),
   this->proposal_parameters[this->covariance_names[i]]);
   }
   */
  return output;
}

double HMCProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                             Particle &old_particle) const
{
  // needs changing
  double output = 0.0;
  /*
   for (size_t i=0; i<this->covariance_names.size(); ++i)
   {
   output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
   old_particle.parameters.get_vector(this->variable_names[i]),
   this->proposal_parameters[this->covariance_names[i]]);
   }
   */
  return output;
}

double HMCProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                         Particle &old_particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  // needs changing
  double output = 0.0;
  /*
  for (size_t i=0; i<this->covariance_names.size(); ++i)
  {
    output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
                              old_particle.parameters.get_vector(this->variable_names[i]),
                              this->proposal_parameters[this->covariance_names[i]]);
  }
  */
  return output;
}

Parameters HMCProposalKernel::simulate(RandomNumberGenerator &rng,
                                       Particle &particle) const
{
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}

Parameters HMCProposalKernel::simulate(RandomNumberGenerator &rng,
                                       Particle &particle,
                                       const Parameters &conditioned_on_parameters) const
{
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 Particle &particle) const
{
  // needs changing
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}
  
Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  // needs changing
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const std::string &variable,
                                                 Particle &particle) const
{
  // needs changing
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const std::string &variable,
                                                 Particle &particle,
                                                 const Parameters &conditioned_on_parameters) const
{
  // needs changing
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  return output;
}

arma::mat HMCProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat HMCProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat HMCProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                     Particle &proposed_particle,
                                                                     Particle &old_particle,
                                                                     const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}
