#include "slice_point_proposal_kernel.h"

namespace ilike
{

SlicePointProposalKernel::SlicePointProposalKernel()
:ProposalKernel(), width(1.0), lower(-arma::datum::inf), upper(arma::datum::inf), max_steps(50)
{
}

SlicePointProposalKernel::~SlicePointProposalKernel()
{
}

SlicePointProposalKernel::SlicePointProposalKernel(const std::string &variable_in, double width_in,
                                                   double lower_in, double upper_in, size_t max_steps_in)
:ProposalKernel(), variable(variable_in), width(width_in), lower(lower_in), upper(upper_in), max_steps(max_steps_in)
{
  if (!(lower_in < upper_in))
    Rcpp::stop("ilike::slice: lower < upper is needed.");
  if (!(width_in > 0.0) && !(std::isfinite(lower_in) && std::isfinite(upper_in)))
    Rcpp::stop("ilike::slice: a positive width is needed unless both lower and upper are finite.");
}

SlicePointProposalKernel::SlicePointProposalKernel(const SlicePointProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void SlicePointProposalKernel::operator=(const SlicePointProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* SlicePointProposalKernel::duplicate() const
{
  return new SlicePointProposalKernel(*this);
}

ProposalKernel* SlicePointProposalKernel::proposal_kernel_duplicate() const
{
  return new SlicePointProposalKernel(*this);
}

void SlicePointProposalKernel::make_copy(const SlicePointProposalKernel &another)
{
  this->variable = another.variable;
  this->width = another.width;
  this->lower = another.lower;
  this->upper = another.upper;
  this->max_steps = another.max_steps;
  this->element = another.element;
  this->value = another.value;
}

double SlicePointProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                         const Particle &old_particle) const
{
  // (not a proposal with a density: SliceMCMC never asks for it)
  Rcpp::stop("ilike::slice is a move of its own (move m = slice(x, ...)), not a proposal.");
}

double SlicePointProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                   const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters SlicePointProposalKernel::simulate(RandomNumberGenerator &rng,
                                                             const Particle &particle) const
{
  // the current point with element `element` of the variable set to `value`
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x = current.get_colvec(this->variable);
  if (this->element >= x.n_elem)
    Rcpp::stop("ilike::slice: element outside the variable.");
  x[this->element] = this->value;
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

Parameters SlicePointProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters SlicePointProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const std::string &variable_in,
                                                                       const Particle &particle) const
{
  if (variable_in != this->variable)
    return Parameters();
  return this->simulate(rng, particle);
}

arma::mat SlicePointProposalKernel::specific_gradient_of_log(const std::string &variable_in,
                                                                            const Particle &proposed_particle,
                                                                            const Particle &old_particle)
{
  Rcpp::stop("ilike::slice: the gradient of the proposal's log density is not available.");
}

arma::mat SlicePointProposalKernel::specific_subsample_gradient_of_log(const std::string &variable_in,
                                                                                      const Particle &proposed_particle,
                                                                                      const Particle &old_particle)
{
  Rcpp::stop("ilike::slice: the gradient of the proposal's log density is not available.");
}

void SlicePointProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* SlicePointProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> SlicePointProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void SlicePointProposalKernel::set_index(Index* index_in)
{
}

void SlicePointProposalKernel::set_index_if_null(Index* index_in)
{
}

bool SlicePointProposalKernel::can_be_evaluated() const
{
  return true;
}

void SlicePointProposalKernel::set_data(Data* data_in)
{
}

}
