#include <iterator>
#include "mirror_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

MirrorProposalKernel::MirrorProposalKernel()
  :SymmetricProposalKernel()
{
}

MirrorProposalKernel::~MirrorProposalKernel()
{
}

// find mean and cov adaptively
MirrorProposalKernel::MirrorProposalKernel(const std::vector<std::string> &variable_names_in)
  :SymmetricProposalKernel()
{
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

// find cov adaptively
MirrorProposalKernel::MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                                           const std::vector<arma::colvec> &means_in)
  :SymmetricProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(means_in[i]);
  }
}

// find mean adaptively
MirrorProposalKernel::MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                                           const std::vector<arma::mat> &covariances_in)
  :SymmetricProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
}

MirrorProposalKernel::MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                                           const std::vector<arma::colvec> &means_in,
                                           const std::vector<arma::mat> &covariances_in)
  :SymmetricProposalKernel()
{
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(means_in[i],covariances_in[i]);
  }
}

MirrorProposalKernel::MirrorProposalKernel(const MirrorProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void MirrorProposalKernel::operator=(const MirrorProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* MirrorProposalKernel::duplicate() const
{
  return( new MirrorProposalKernel(*this));
}

ProposalKernel* MirrorProposalKernel::proposal_kernel_duplicate() const
{
  return( new MirrorProposalKernel(*this));
}

SymmetricProposalKernel* MirrorProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new MirrorProposalKernel(*this));
}

void MirrorProposalKernel::make_copy(const MirrorProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
}

double MirrorProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_vector(i->first),
                                            2.0*mean-old_particle.move_parameters->get_vector(i->first),
                                            (1.0/sqrt(scale))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
  
  // meed to write something a bit like get_grad, since need a flag to say if mcmc_trans_params is already set
  // need to have piece of code that does the transform before prop, then does it back again, also gathers Jacobian - in ProposalKernel since its the same piece of code everywhere. Maybe also has a pointer to switch params to transformed space, then proposal does not need if statement in simulate. What about gradient bit though? Needs to know what space we are in. Particle can deal with this as well, since it is storing the info about whether we are in a transfoemed space.
  
}

double MirrorProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

double MirrorProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

double MirrorProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

arma::mat MirrorProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

Parameters MirrorProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle) const
{
  Parameters output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    output[i->first] = rmvnorm(rng,
                               2.0*i->second.get_mean()-particle.move_parameters->get_vector(i->first),
                               sqrt(i->second.get_double_scale())*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters MirrorProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng,particle);
}

Parameters MirrorProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters MirrorProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters MirrorProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output = *particle.move_parameters;
  output[variable] = rmvnorm(rng,
                             2.0*found->second.get_mean()-particle.move_parameters->get_vector(variable),
                             sqrt(found->second.get_double_scale())*found->second.get_chol(),
                             true);
  return output;
}

Parameters MirrorProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->subsample_simulate(rng,
                                  variable,
                                  particle);
}

arma::mat MirrorProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                     Particle &proposed_particle,
                                                                     Particle &old_particle)
{
  Rcpp::stop("MirrorProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat MirrorProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                     Particle &proposed_particle,
                                                                     Particle &old_particle,
                                                                     const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("MirrorProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat MirrorProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                               Particle &proposed_particle,
                                                                               Particle &old_particle,
                                                                               const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("MirrorProposalKernel::specific_gradient_of_log - not written yet.");
}
