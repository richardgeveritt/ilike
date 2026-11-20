#include "choice_proposal_kernel.h"

namespace ilike
{

// --- ChoiceProposalKernel ----------------------------------------------------

ChoiceProposalKernel::ChoiceProposalKernel()
:ProposalKernel()
{
}

ChoiceProposalKernel::ChoiceProposalKernel(const std::vector<ProposalKernel*> &components_in)
:ProposalKernel()
{
  for (ProposalKernel* component : components_in)
    this->components.push_back(std::unique_ptr<ProposalKernel>(component));
}

ChoiceProposalKernel::~ChoiceProposalKernel()
{
}

ChoiceProposalKernel::ChoiceProposalKernel(const ChoiceProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void ChoiceProposalKernel::operator=(const ChoiceProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* ChoiceProposalKernel::duplicate() const
{
  return new ChoiceProposalKernel(*this);
}

ProposalKernel* ChoiceProposalKernel::proposal_kernel_duplicate() const
{
  return new ChoiceProposalKernel(*this);
}

void ChoiceProposalKernel::make_copy(const ChoiceProposalKernel &another)
{
  this->components.clear();
  for (const auto &component : another.components)
    this->components.push_back(std::unique_ptr<ProposalKernel>(component->proposal_kernel_duplicate()));
}

const ProposalKernel* ChoiceProposalKernel::chosen(const Particle &particle) const
{
  const double* value = this->tuning_values(particle);
  if (value == NULL)
    Rcpp::stop("choose([...]) is used in SMC_with_MCMC, where each particle carries the kernel it uses.");
  long k = std::lround(value[0]);
  if ((k < 0) || (size_t(k) >= this->components.size()))
    Rcpp::stop("choose: a particle's kernel index is out of range.");
  return this->components[size_t(k)].get();
}

Particle ChoiceProposalKernel::move(RandomNumberGenerator &rng,
                                    const Particle &particle,
                                    const Index* new_outputs_index) const
{
  return this->chosen(particle)->move(rng, particle, new_outputs_index);
}

void ChoiceProposalKernel::move_into(RandomNumberGenerator &rng,
                                     const Particle &particle,
                                     const Index* new_outputs_index,
                                     Particle &scratch) const
{
  this->chosen(particle)->move_into(rng, particle, new_outputs_index, scratch);
}

Parameters ChoiceProposalKernel::simulate(RandomNumberGenerator &rng,
                                          const Particle &particle) const
{
  return this->chosen(particle)->simulate(rng, particle);
}

Parameters ChoiceProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const Particle &particle) const
{
  Rcpp::stop("choose with subsampling is not available.");
}

Parameters ChoiceProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    const Particle &particle) const
{
  Rcpp::stop("choose with subsampling is not available.");
}

double ChoiceProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                      const Particle &old_particle) const
{
  // (the two particles carry the same tuning, so the same component)
  return this->chosen(old_particle)->evaluate_kernel(proposed_particle, old_particle);
}

double ChoiceProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                const Particle &old_particle) const
{
  Rcpp::stop("choose with subsampling is not available.");
}

arma::mat ChoiceProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         const Particle &proposed_particle,
                                                         const Particle &old_particle)
{
  Rcpp::stop("choose: the gradient of the proposal's log density is not available.");
}

arma::mat ChoiceProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                   const Particle &proposed_particle,
                                                                   const Particle &old_particle)
{
  Rcpp::stop("choose: the gradient of the proposal's log density is not available.");
}

void ChoiceProposalKernel::smc_adapt(SMCOutput* current_state)
{
  ProposalKernel::smc_adapt(current_state);
  for (auto &component : this->components)
    component->smc_adapt(current_state);
}

void ChoiceProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto &component : this->components)
    component->set_proposal_parameters(proposal_parameters_in);
}

GradientEstimatorOutput* ChoiceProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> ChoiceProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output(1, this);
  for (const auto &component : this->components)
  {
    std::vector<const ProposalKernel*> inside = component->get_proposals();
    output.insert(output.end(), inside.begin(), inside.end());
  }
  return output;
}

std::vector<const ProposalKernel*> ChoiceProposalKernel::proposals_used(const Particle &particle) const
{
  std::vector<const ProposalKernel*> output(1, this);
  const double* value = this->tuning_values(particle);
  if (value == NULL)
    return output;
  long k = std::lround(value[0]);
  if ((k >= 0) && (size_t(k) < this->components.size()))
  {
    std::vector<const ProposalKernel*> inside = this->components[size_t(k)]->get_proposals();
    output.insert(output.end(), inside.begin(), inside.end());
  }
  return output;
}

void ChoiceProposalKernel::set_index(Index* index_in)
{
  for (auto &component : this->components)
    component->set_index(index_in);
}

void ChoiceProposalKernel::set_index_if_null(Index* index_in)
{
  for (auto &component : this->components)
    component->set_index_if_null(index_in);
}

bool ChoiceProposalKernel::can_be_evaluated() const
{
  for (const auto &component : this->components)
    if (!component->can_be_evaluated())
      return false;
  return true;
}

void ChoiceProposalKernel::set_data(Data* data_in)
{
  for (auto &component : this->components)
    component->set_data(data_in);
}

// --- ChoiceAdaptor -----------------------------------------------------------

ChoiceAdaptor::ChoiceAdaptor()
:ESJDAdaptor(), change(0.1)
{
}

ChoiceAdaptor::ChoiceAdaptor(const std::string &variable_in, const std::vector<double> &costs_in, double change_in,
                             bool full_metric_in)
:ESJDAdaptor(variable_in, {"kernel"}, {0.0}, {double(costs_in.size()) - 1.0}, {true}, 0.0, 0.0, full_metric_in, -1),
 costs(costs_in), change(change_in)
{
}

ChoiceAdaptor::~ChoiceAdaptor()
{
}

ChoiceAdaptor::ChoiceAdaptor(const ChoiceAdaptor &another)
:ESJDAdaptor(another)
{
  this->make_copy(another);
}

void ChoiceAdaptor::operator=(const ChoiceAdaptor &another)
{
  if (this == &another)
    return;
  ESJDAdaptor::operator=(another);
  this->make_copy(another);
}

SMCAdaptor* ChoiceAdaptor::duplicate() const
{
  return new ChoiceAdaptor(*this);
}

void ChoiceAdaptor::make_copy(const ChoiceAdaptor &another)
{
  this->costs = another.costs;
  this->change = another.change;
}

double ChoiceAdaptor::perturbed(RandomNumberGenerator &rng, size_t k, double v) const
{
  size_t K = this->costs.size();
  if ((K < 2) || !(runif(rng) < this->change))
    return v;
  // another kernel, uniformly
  size_t other = std::min(size_t(std::floor(runif(rng)*double(K - 1))), K - 2);
  size_t current = size_t(std::lround(v));
  return double((other >= current) ? other + 1 : other);
}

void ChoiceAdaptor::score(size_t i, const Parameters &before, const Parameters &after, const ParticleTuning &tuning)
{
  if (i >= this->scores.n_elem)
    return;
  const ParticleTuning::Row* row = NULL;
  for (const auto &r : tuning.rows)
    if (r.adaptor == this)
      row = &r;
  double distance = ((row != NULL) && (row->iterations > 0)) ? row->jumps/double(row->iterations)
                                                             : this->squared_jump(before, after);
  size_t k = size_t(std::lround(this->values(0, i)));
  double cost = (k < this->costs.size()) ? this->costs[k] : 1.0;
  distance = distance/std::max(1e-12, cost);
  this->scores[i] = std::isfinite(distance) ? distance : 0.0;
}

void ChoiceAdaptor::summary(std::vector<std::string> &names_out, std::vector<double> &numbers) const
{
  // the share of the particles using each kernel
  for (size_t k = 0; k < this->costs.size(); ++k)
  {
    double share = 0.0;
    for (arma::uword i = 0; i < this->values.n_cols; ++i)
      if (size_t(std::lround(this->values(0, i))) == k)
        share = share + 1.0;
    names_out.push_back("kernel" + std::to_string(k + 1) + ".share");
    numbers.push_back((this->values.n_cols > 0) ? share/double(this->values.n_cols) : arma::datum::nan);
  }
}

}
