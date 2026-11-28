#include "involution_proposal_kernel.h"

namespace ilike
{

InvolutionProposalKernel::InvolutionProposalKernel()
  :ProposalKernel()
{
}

InvolutionProposalKernel::~InvolutionProposalKernel()
{
}

InvolutionProposalKernel::InvolutionProposalKernel(const InvolutionProposalKernel &another)
  :ProposalKernel(another)
{
}

void InvolutionProposalKernel::operator=(const InvolutionProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
}

Kernel* InvolutionProposalKernel::duplicate() const
{
  return new InvolutionProposalKernel(*this);
}

ProposalKernel* InvolutionProposalKernel::proposal_kernel_duplicate() const
{
  return new InvolutionProposalKernel(*this);
}

double InvolutionProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                          const Particle &old_particle) const
{
  return 0.0;
}

double InvolutionProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                    const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters InvolutionProposalKernel::simulate(RandomNumberGenerator &rng,
                                              const Particle &particle) const
{
  return Parameters();
}

Parameters InvolutionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                        const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters InvolutionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                        const std::string &variable,
                                                        const Particle &particle) const
{
  return this->simulate(rng, particle);
}

arma::mat InvolutionProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                             const Particle &proposed_particle,
                                                             const Particle &old_particle)
{
  Rcpp::stop("involution(…): the gradient of a proposal's log density is not available.");
}

arma::mat InvolutionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                       const Particle &proposed_particle,
                                                                       const Particle &old_particle)
{
  Rcpp::stop("involution(…): the gradient of a proposal's log density is not available.");
}

void InvolutionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* InvolutionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> InvolutionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void InvolutionProposalKernel::set_index(Index* index_in)
{
}

void InvolutionProposalKernel::set_index_if_null(Index* index_in)
{
}

bool InvolutionProposalKernel::can_be_evaluated() const
{
  return true;
}

void InvolutionProposalKernel::set_data(Data* data_in)
{
}

}
