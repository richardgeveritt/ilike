#include "esjd_adaptor.h"
#include "covariance_tools.h"
#include "proposal_kernel.h"
#include "smc_output.h"
#include "transform.h"
#include "particle.h"

namespace ilike
{

// --- ESJDAdaptor -------------------------------------------------------------

ESJDAdaptor::ESJDAdaptor()
:SMCAdaptor(), constant(0.0), perturb(0.05), full_metric(false), cost_index(-1)
{
}

ESJDAdaptor::ESJDAdaptor(const std::string &variable_in,
                         const std::vector<std::string> &names_in,
                         const std::vector<double> &lower_in,
                         const std::vector<double> &upper_in,
                         const std::vector<bool> &integer_in,
                         double constant_in,
                         double perturb_in,
                         bool full_metric_in,
                         int cost_index_in)
:SMCAdaptor(), variable(variable_in), names(names_in), lower(lower_in), upper(upper_in), integer(integer_in),
 constant(constant_in), perturb(perturb_in), full_metric(full_metric_in), cost_index(cost_index_in)
{
}

ESJDAdaptor::~ESJDAdaptor()
{
}

ESJDAdaptor::ESJDAdaptor(const ESJDAdaptor &another)
:SMCAdaptor(another)
{
  this->make_copy(another);
}

void ESJDAdaptor::operator=(const ESJDAdaptor &another)
{
  if (this == &another)
    return;
  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

SMCAdaptor* ESJDAdaptor::duplicate() const
{
  return new ESJDAdaptor(*this);
}

void ESJDAdaptor::make_copy(const ESJDAdaptor &another)
{
  this->variable = another.variable;
  this->names = another.names;
  this->lower = another.lower;
  this->upper = another.upper;
  this->integer = another.integer;
  this->constant = another.constant;
  this->perturb = another.perturb;
  this->full_metric = another.full_metric;
  this->cost_index = another.cost_index;
  this->values = another.values;
  this->scores = another.scores;
  this->metric = another.metric;
  this->updates = another.updates;
}

void ESJDAdaptor::smc_adapt(SMCOutput* current_state)
{
  // (the work is done by update, which the SMC move calls with its stream)
}

void ESJDAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  Rcpp::stop("smc_adapt = esjd(…) is for SMC_with_MCMC, not for ensemble Kalman methods.");
}

void ESJDAdaptor::initialise(RandomNumberGenerator &rng, size_t number_of_particles)
{
  size_t p = this->names.size();
  this->values.set_size(p, number_of_particles);
  for (size_t i = 0; i < number_of_particles; ++i)
    for (size_t k = 0; k < p; ++k)
    {
      double u = runif(rng);
      if (this->integer[k])
        this->values(k, i) = std::floor(this->lower[k] + u*(this->upper[k] - this->lower[k] + 1.0));
      else
        this->values(k, i) = this->lower[k] + u*(this->upper[k] - this->lower[k]);
    }
}

void ESJDAdaptor::update(RandomNumberGenerator &rng, SMCOutput* current_state, size_t number_of_particles)
{
  size_t N = number_of_particles;
  size_t p = this->names.size();
  this->updates = this->updates + 1;
  // (particles not scored, NaN, have weight 0: in choose, those whose
  // particle used another kernel, C224)
  arma::uvec finite = arma::find_finite(this->scores);
  bool scored = (this->values.n_cols == N) && (this->scores.n_elem == N) && (finite.n_elem > 0);
  if (this->values.n_cols != N)
    this->initialise(rng, N);
  else if (scored)
  {
    // resample with weights a + Lambda (systematic), perturb, shuffle
    arma::colvec weights(N, arma::fill::zeros);
    weights.elem(finite) = this->constant + arma::clamp(arma::colvec(this->scores.elem(finite)), 0.0, arma::datum::inf);
    double total = arma::accu(weights);
    if (!(total > 0.0) || !std::isfinite(total))
      weights.fill(1.0/double(N));
    else
      weights = weights/total;
    arma::colvec cumulative = arma::cumsum(weights);
    std::vector<size_t> chosen(N);
    double u = runif(rng);
    size_t j = 0;
    for (size_t i = 0; i < N; ++i)
    {
      double point = (double(i) + u)/double(N);
      while ((j + 1 < N) && (cumulative[j] < point))
        ++j;
      chosen[i] = j;
    }
    for (size_t i = N; i > 1; --i)
    {
      size_t k = size_t(std::floor(runif(rng)*double(i)));
      if (k >= i)
        k = i - 1;
      std::swap(chosen[i - 1], chosen[k]);
    }
    arma::mat next(p, N);
    for (size_t i = 0; i < N; ++i)
      for (size_t k = 0; k < p; ++k)
        next(k, i) = this->perturbed(rng, k, this->values(k, chosen[i]));
    this->values = next;
  }
  this->scores.set_size(N);
  this->scores.fill(arma::datum::nan);
  this->compute_metric(current_state);
}

double ESJDAdaptor::perturbed(RandomNumberGenerator &rng, size_t k, double v) const
{
  if (this->lower[k] == this->upper[k])
    return this->lower[k]; // (a fixed value: not tuned)
  if (this->integer[k])
  {
    double step = std::floor(runif(rng)*3.0) - 1.0;
    return std::max(1.0, v + std::min(1.0, step));
  }
  return v*std::exp(this->perturb*rnorm(rng));
}

void ESJDAdaptor::compute_metric(SMCOutput* current_state)
{
  // the metric: the particles' variances (or covariance) of the block, on the
  // kernel's scale
  Particles &generation = current_state->back();
  std::vector<std::string> block = split_stacked_name(this->variable);
  std::vector<size_t> offsets;
  const Transform* transform = (this->proposal != NULL) ? this->proposal->get_transform() : NULL;
  arma::mat points = population_points(generation, block, offsets, transform);
  arma::colvec weights;
  size_t K = split_stacked_name(this->variable).size();
  size_t d = (offsets.size() > K) ? offsets[K] : points.n_cols;
  if ((points.n_rows < 2) || !population_weights(generation, points.n_rows, weights))
  {
    this->metric = arma::ones<arma::colvec>(d);
    this->full_metric = false;
    return;
  }
  arma::mat x = points.cols(0, d - 1);
  arma::rowvec mean = weights.t()*x;
  arma::mat centred = x.each_row() - mean;
  if (this->full_metric)
  {
    arma::mat covariance = centred.t()*(centred.each_col() % weights);
    arma::mat inverse;
    if (arma::inv_sympd(inverse, arma::symmatu(covariance)))
    {
      this->metric = inverse;
      return;
    }
  }
  arma::colvec variances = arma::vectorise((centred % centred).t()*weights);
  variances.transform([](double v) { return (v > 0.0) ? 1.0/v : 1.0; });
  this->metric = variances;
}

void ESJDAdaptor::score(size_t i, const Parameters &before, const Parameters &after, const ParticleTuning &tuning)
{
  if (i >= this->scores.n_elem)
    return;
  const ParticleTuning::Row* row = NULL;
  for (const auto &r : tuning.rows)
    if (r.adaptor == this)
      row = &r;
  // in choose (C224), a kernel the particle did not use is not scored
  if ((row != NULL) && (row->iterations == 0) && tuning.has_choice())
  {
    this->scores[i] = arma::datum::nan;
    return;
  }
  double distance = ((row != NULL) && (row->iterations > 0)) ? row->jumps/double(row->iterations)
                                                             : this->squared_jump(before, after);
  if (this->cost_index >= 0)
    distance = distance/std::max(1.0, this->values(size_t(this->cost_index), i));
  this->scores[i] = std::isfinite(distance) ? distance : 0.0;
}

void observe_tuned_proposal(const Particle &current,
                            const Particle &proposed,
                            const ProposalKernel* proposal,
                            double alpha)
{
  const ParticleTuning* tuning = current.tuning;
  if ((tuning == NULL) || (proposal == NULL))
    return;
  std::vector<const ProposalKernel*> kernels = proposal->proposals_used(current);
  for (const auto &row : tuning->rows)
  {
    bool within = false;
    for (const ProposalKernel* kernel : kernels)
      if ((kernel != NULL) && (kernel->get_instance_index() == row.instance))
        within = true;
    if (!within)
      continue;
    double jump = (alpha > 0.0) ? alpha*row.adaptor->squared_jump(current.parameters, proposed.parameters) : 0.0;
    row.jumps = row.jumps + (std::isfinite(jump) ? jump : 0.0);
    row.iterations = row.iterations + 1;
  }
}

double ESJDAdaptor::squared_jump(const Parameters &before, const Parameters &after) const
{
  const Transform* transform = (this->proposal != NULL) ? this->proposal->get_transform() : NULL;
  arma::colvec x0 = (transform != NULL) ? transform->transform(before).get_colvec(this->variable) : before.get_colvec(this->variable);
  arma::colvec x1 = (transform != NULL) ? transform->transform(after).get_colvec(this->variable) : after.get_colvec(this->variable);
  return this->mahalanobis(x1 - x0);
}

double ESJDAdaptor::mahalanobis(const arma::colvec &jump) const
{
  if ((this->metric.n_cols == jump.n_elem) && (this->metric.n_rows == jump.n_elem))
    return arma::as_scalar(jump.t()*this->metric*jump);
  if (this->metric.n_elem == jump.n_elem)
    return arma::accu(jump % jump % arma::vectorise(this->metric));
  return arma::dot(jump, jump);
}

void ESJDAdaptor::summary(std::vector<std::string> &names_out, std::vector<double> &numbers) const
{
  static const char* statistics[] = {"mean", "sd", "q10", "q50", "q90"};
  for (size_t k = 0; k < this->names.size(); ++k)
  {
    arma::rowvec row = (this->values.n_cols > 0) ? arma::rowvec(this->values.row(k)) : arma::rowvec();
    arma::colvec q = (row.n_elem > 0) ? arma::colvec(arma::quantile(row, arma::vec({0.1, 0.5, 0.9})).t())
                                      : arma::colvec(3, arma::fill::value(arma::datum::nan));
    double m = (row.n_elem > 0) ? arma::mean(row) : arma::datum::nan;
    double s = (row.n_elem > 1) ? arma::stddev(row) : arma::datum::nan;
    double values_out[] = {m, s, q[0], q[1], q[2]};
    for (size_t j = 0; j < 5; ++j)
    {
      names_out.push_back(this->names[k] + "." + statistics[j]);
      numbers.push_back(values_out[j]);
    }
  }
}

// --- CompositeSMCAdaptor -----------------------------------------------------

CompositeSMCAdaptor::CompositeSMCAdaptor()
:SMCAdaptor()
{
}

CompositeSMCAdaptor::CompositeSMCAdaptor(const std::vector<SMCAdaptor*> &parts_in)
:SMCAdaptor()
{
  for (SMCAdaptor* part : parts_in)
    this->parts.push_back(std::unique_ptr<SMCAdaptor>(part));
}

CompositeSMCAdaptor::~CompositeSMCAdaptor()
{
}

CompositeSMCAdaptor::CompositeSMCAdaptor(const CompositeSMCAdaptor &another)
:SMCAdaptor(another)
{
  this->make_copy(another);
}

void CompositeSMCAdaptor::operator=(const CompositeSMCAdaptor &another)
{
  if (this == &another)
    return;
  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

SMCAdaptor* CompositeSMCAdaptor::duplicate() const
{
  return new CompositeSMCAdaptor(*this);
}

void CompositeSMCAdaptor::make_copy(const CompositeSMCAdaptor &another)
{
  this->parts.clear();
  for (const auto &part : another.parts)
    this->parts.push_back(std::unique_ptr<SMCAdaptor>(part->duplicate()));
}

void CompositeSMCAdaptor::set_proposal(ProposalKernel* proposal_in)
{
  SMCAdaptor::set_proposal(proposal_in);
  for (auto &part : this->parts)
    part->set_proposal(proposal_in);
}

void CompositeSMCAdaptor::smc_adapt(SMCOutput* current_state)
{
  for (auto &part : this->parts)
    part->smc_adapt(current_state);
}

void CompositeSMCAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  for (auto &part : this->parts)
    part->ensemble_adapt(current_state);
}

std::vector<SMCAdaptor*> CompositeSMCAdaptor::get_parts() const
{
  std::vector<SMCAdaptor*> output;
  for (const auto &part : this->parts)
    output.push_back(part.get());
  return output;
}

}
