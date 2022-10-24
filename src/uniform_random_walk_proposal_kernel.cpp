#include <iterator>
#include "uniform_random_walk_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel()
  :SymmetricProposalKernel()
{
  this->unused_variables_kept = true;
}

UniformRandomWalkProposalKernel::~UniformRandomWalkProposalKernel()
{
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in)
  :SymmetricProposalKernel()
{
  this->unused_variables_kept = true;
  
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = arma::mat();
  }
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                 const std::vector<arma::mat> &halfwidths_in)
  :SymmetricProposalKernel()
{
  this->unused_variables_kept = true;
  
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = halfwidths_in[i];
  }
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const UniformRandomWalkProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void UniformRandomWalkProposalKernel::operator=(const UniformRandomWalkProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* UniformRandomWalkProposalKernel::duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

ProposalKernel* UniformRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

SymmetricProposalKernel* UniformRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

void UniformRandomWalkProposalKernel::make_copy(const UniformRandomWalkProposalKernel &another)
{
  this->unused_variables_kept = another.unused_variables_kept;
  this->proposal_info = another.proposal_info;
}

double UniformRandomWalkProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                                 Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    //output = output - double(i->second.n_elem)*log(2.0) - arma::accu(arma::log(i->second));
    arma::mat current_position = (*old_particle.move_parameters)[i->first];
    output = output + dunif((*proposed_particle.move_parameters)[i->first],
                            current_position-i->second,
                            current_position+i->second);
  }
  return output;
}

double UniformRandomWalkProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                                 Particle &old_particle,
                                                                 const Parameters &conditioned_on_parameters) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

double UniformRandomWalkProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                           Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

double UniformRandomWalkProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                         Particle &old_particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters UniformRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                      Particle &particle) const
{
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::mat current_position = (*particle.move_parameters)[i->first];
    output[i->first] = runif(rng,
                             current_position-i->second,
                             current_position+i->second);
  }
  return output;
}

Parameters UniformRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                     Particle &particle,
                                                     const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               Particle &particle,
                                                               const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const std::string &variable,
                                                               Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found =  this->proposal_info.find(variable);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  arma::mat current_position = (*particle.move_parameters)[found->first];
  output[found->first] = runif(rng,
                               current_position-found->second,
                               current_position+found->second);
  return output;
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const std::string &variable,
                                                                Particle &particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->subsample_simulate(rng,
                                  variable,
                                  particle);
}

arma::mat UniformRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                    Particle &proposed_particle,
                                                                    Particle &old_particle)
{
  Rcpp::stop("UniformRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat UniformRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                    Particle &proposed_particle,
                                                                    Particle &old_particle,
                                                                    const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("UniformRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat UniformRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                              Particle &proposed_particle,
                                                                              Particle &old_particle,
                                                                              const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("UniformRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}
