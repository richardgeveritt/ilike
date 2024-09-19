#include <iterator>
#include "gamma_independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

namespace ilike
{
GammaIndependentProposalKernel::GammaIndependentProposalKernel()
  :IndependentProposalKernel()
{
}

GammaIndependentProposalKernel::~GammaIndependentProposalKernel()
{
}

GammaIndependentProposalKernel::GammaIndependentProposalKernel(const std::string &variable_in,
                                                               const double &shape_in,
                                                               const double &rate_in)
:IndependentProposalKernel()
{
  this->variable = variable_in;
  this->shape = shape_in;
  this->rate = rate_in;
}

GammaIndependentProposalKernel::GammaIndependentProposalKernel(const GammaIndependentProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void GammaIndependentProposalKernel::operator=(const GammaIndependentProposalKernel &another)
{
  if(this == &another)
    return;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* GammaIndependentProposalKernel::duplicate() const
{
  return( new GammaIndependentProposalKernel(*this));
}

ProposalKernel* GammaIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new GammaIndependentProposalKernel(*this));
}

IndependentProposalKernel* GammaIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new GammaIndependentProposalKernel(*this));
}

void GammaIndependentProposalKernel::make_copy(const GammaIndependentProposalKernel &another)
{
  this->variable = another.variable;
  this->shape = another.shape;
  this->rate = another.rate;
}

double GammaIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  return dgamma(proposed_particle[this->variable][0],
                this->shape,
                this->rate);
}

/*
double GammaIndependentProposalKernel::evaluate_independent_kernel(Variables* proposed_particle,
                                                                      const Parameters &conditioned_on_parameters) const
{
  return this->evaluate_independent_kernel(proposed_particle);
}
*/

double GammaIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->evaluate_independent_kernel(proposed_particle);
}

Parameters GammaIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  Parameters output;
  output[this->variable] = rgamma(rng, this->shape, this->rate);
  return output;
}

Parameters GammaIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng);
}

Parameters GammaIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters GammaIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                          const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->independent_simulate(rng);
}

Parameters GammaIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                          const std::string &variable_in) const
{
  Parameters output;
  if (this->variable==variable_in)
    output[this->variable] = rgamma(rng, this->shape, this->rate);
  return output;
}

Parameters GammaIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                          const std::string &variable_in,
                                                                          const Parameters &conditioned_on_parameters) const
{
  Parameters output;
  if (this->variable==variable_in)
    output[this->variable] = rgamma(rng, this->shape, this->rate);
  return output;
}

arma::mat GammaIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                      const Parameters &proposed_particle)
{
  Rcpp::stop("GammaIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

arma::mat GammaIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                const Parameters &proposed_particle)
{
  Rcpp::stop("GammaIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

void GammaIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

GradientEstimatorOutput* GammaIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> GammaIndependentProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void GammaIndependentProposalKernel::set_index(Index* index_in)
{
}

void GammaIndependentProposalKernel::set_index_if_null(Index* index_in)
{
}

bool GammaIndependentProposalKernel::can_be_evaluated() const
{
  return true;
}

void GammaIndependentProposalKernel::set_data(Data* data_in)
{
  
}
}
