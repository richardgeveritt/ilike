#include "dimension_jump_proposal_kernel.h"

namespace ilike
{

DimensionJumpProposalKernel::DimensionJumpProposalKernel()
:ProposalKernel(), minimum(1.0), maximum(1.0)
{
}

DimensionJumpProposalKernel::~DimensionJumpProposalKernel()
{
}

DimensionJumpProposalKernel::DimensionJumpProposalKernel(const std::string &variable_in, double minimum_in, double maximum_in)
:ProposalKernel(), variable(variable_in), minimum(minimum_in), maximum(maximum_in)
{
  if (!(minimum_in < maximum_in))
    Rcpp::stop("dimension_jump: " + variable_in + " needs at least two values (minimum < maximum).");
}

DimensionJumpProposalKernel::DimensionJumpProposalKernel(const DimensionJumpProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void DimensionJumpProposalKernel::operator=(const DimensionJumpProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* DimensionJumpProposalKernel::duplicate() const
{
  return new DimensionJumpProposalKernel(*this);
}

ProposalKernel* DimensionJumpProposalKernel::proposal_kernel_duplicate() const
{
  return new DimensionJumpProposalKernel(*this);
}

void DimensionJumpProposalKernel::make_copy(const DimensionJumpProposalKernel &another)
{
  this->variable = another.variable;
  this->minimum = another.minimum;
  this->maximum = another.maximum;
}

double DimensionJumpProposalKernel::log_density(double to, double from) const
{
  if (std::abs(std::abs(to - from) - 1.0) > 1e-9)
    return -arma::datum::inf;
  if ((to < this->minimum) || (to > this->maximum))
    return -arma::datum::inf;
  if ((from <= this->minimum) || (from >= this->maximum))
    return 0.0; // (forced at a bound)
  return std::log(0.5);
}

double DimensionJumpProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                             const Particle &old_particle) const
{
  return this->log_density(proposed_particle.parameters[this->variable](0), old_particle.parameters[this->variable](0));
}

double DimensionJumpProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                       const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters DimensionJumpProposalKernel::simulate(RandomNumberGenerator &rng,
                                                 const Particle &particle) const
{
  const double k = particle.parameters[this->variable](0);
  double to;
  if (k <= this->minimum)
    to = k + 1.0;
  else if (k >= this->maximum)
    to = k - 1.0;
  else
    to = (runif(rng) < 0.5) ? k + 1.0 : k - 1.0;
  Parameters output;
  output[this->variable] = arma::mat(1, 1, arma::fill::value(to));
  return output;
}

void DimensionJumpProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                const Particle &particle,
                                                ProposalTarget &target) const
{
  target.write_all(this->simulate(rng, particle));
}

Parameters DimensionJumpProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                           const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters DimensionJumpProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                           const std::string &variable_in,
                                                           const Particle &particle) const
{
  if (variable_in != this->variable)
    return Parameters();
  return this->simulate(rng, particle);
}

arma::mat DimensionJumpProposalKernel::specific_gradient_of_log(const std::string &variable_in,
                                                                const Particle &proposed_particle,
                                                                const Particle &old_particle)
{
  Rcpp::stop("dimension_jump: the proposal is discrete; it has no gradient.");
}

arma::mat DimensionJumpProposalKernel::specific_subsample_gradient_of_log(const std::string &variable_in,
                                                                          const Particle &proposed_particle,
                                                                          const Particle &old_particle)
{
  Rcpp::stop("dimension_jump: the proposal is discrete; it has no gradient.");
}

void DimensionJumpProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* DimensionJumpProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> DimensionJumpProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void DimensionJumpProposalKernel::set_index(Index* index_in)
{
}

void DimensionJumpProposalKernel::set_index_if_null(Index* index_in)
{
}

bool DimensionJumpProposalKernel::can_be_evaluated() const
{
  return true;
}

void DimensionJumpProposalKernel::set_data(Data* data_in)
{
}

}
