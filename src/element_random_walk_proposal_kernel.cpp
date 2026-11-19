#include "element_random_walk_proposal_kernel.h"
#include "proposal_target.h"

namespace ilike
{

ElementRandomWalkProposalKernel::ElementRandomWalkProposalKernel()
  :SymmetricProposalKernel()
{
  this->row = 0;
  this->column = 0;
  this->scale = 1.0;
  this->step_kind = ElementStep::normal;
}

ElementRandomWalkProposalKernel::ElementRandomWalkProposalKernel(const std::string &variable_in,
                                                                 size_t row_in,
                                                                 size_t column_in,
                                                                 double scale_in,
                                                                 ElementStep step_in)
  :SymmetricProposalKernel()
{
  this->variable = variable_in;
  this->row = row_in;
  this->column = column_in;
  this->scale = scale_in;
  this->step_kind = step_in;
}

ElementRandomWalkProposalKernel::~ElementRandomWalkProposalKernel()
{
}

ElementRandomWalkProposalKernel::ElementRandomWalkProposalKernel(const ElementRandomWalkProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void ElementRandomWalkProposalKernel::operator=(const ElementRandomWalkProposalKernel &another)
{
  if (this == &another)
    return;
  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* ElementRandomWalkProposalKernel::duplicate() const
{
  return new ElementRandomWalkProposalKernel(*this);
}

ProposalKernel* ElementRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return new ElementRandomWalkProposalKernel(*this);
}

SymmetricProposalKernel* ElementRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return new ElementRandomWalkProposalKernel(*this);
}

void ElementRandomWalkProposalKernel::make_copy(const ElementRandomWalkProposalKernel &another)
{
  this->variable = another.variable;
  this->row = another.row;
  this->column = another.column;
  this->scale = another.scale;
  this->step_kind = another.step_kind;
}

const arma::mat &ElementRandomWalkProposalKernel::current_value(const Particle &particle) const
{
  const arma::mat &value = (particle.get_transformed_parameters(this))[this->variable];
  if ((this->row >= value.n_rows) || (this->column >= value.n_cols))
    Rcpp::stop("The move on " + this->variable + "[" + std::to_string(this->row + 1) +
               (value.n_cols > 1 || this->column > 0 ? ", " + std::to_string(this->column + 1) : std::string("")) +
               "] is outside the variable, which has " + std::to_string(value.n_rows) + " row(s) and " +
               std::to_string(value.n_cols) + " column(s).");
  return value;
}

double ElementRandomWalkProposalKernel::step(RandomNumberGenerator &rng, double current) const
{
  if (this->step_kind == ElementStep::normal)
    return rnorm(rng, current, this->scale);
  if (this->step_kind == ElementStep::uniform)
    return runif(rng, current - this->scale, current + this->scale);
  return 1.0 - current;
}

double ElementRandomWalkProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                 const Particle &old_particle) const
{
  double from = this->current_value(old_particle)(this->row, this->column);
  double to = this->current_value(proposed_particle)(this->row, this->column);
  if (this->step_kind == ElementStep::normal)
    return dnorm(to, from, this->scale);
  if (this->step_kind == ElementStep::uniform)
    return dunif(to, from - this->scale, from + this->scale);
  return 0.0; // (flip: the reverse move is certain too)
}

double ElementRandomWalkProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                           const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

void ElementRandomWalkProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                    const Particle &particle,
                                                    ProposalTarget &target) const
{
  const arma::mat &current = this->current_value(particle);
  double proposed = this->step(rng, current(this->row, this->column));
  arma::mat &written = target.matrix(this->variable);
  written = current;
  written(this->row, this->column) = proposed;
}

Parameters ElementRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                     const Particle &particle) const
{
  const arma::mat &current = this->current_value(particle);
  double proposed = this->step(rng, current(this->row, this->column));
  Parameters output;
  output[this->variable] = current;
  output[this->variable](this->row, this->column) = proposed;
  return output;
}

Parameters ElementRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                     const Particle &particle,
                                                     const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng, particle);
}

Parameters ElementRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters ElementRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const Particle &particle,
                                                               const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng, particle);
}

Parameters ElementRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const std::string &variable,
                                                               const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters ElementRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const std::string &variable,
                                                               const Particle &particle,
                                                               const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng, particle);
}

arma::mat ElementRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                    const Particle &proposed_particle,
                                                                    const Particle &old_particle)
{
  Rcpp::stop("ElementRandomWalkProposalKernel: no gradient.");
}

arma::mat ElementRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                              const Particle &proposed_particle,
                                                                              const Particle &old_particle)
{
  Rcpp::stop("ElementRandomWalkProposalKernel: no gradient.");
}

void ElementRandomWalkProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* ElementRandomWalkProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> ElementRandomWalkProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void ElementRandomWalkProposalKernel::set_index(Index* index_in)
{
}

void ElementRandomWalkProposalKernel::set_index_if_null(Index* index_in)
{
}

bool ElementRandomWalkProposalKernel::can_be_evaluated() const
{
  return true;
}

void ElementRandomWalkProposalKernel::set_data(Data* data_in)
{
}

}
