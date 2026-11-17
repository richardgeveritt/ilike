#include "plate_likelihood_estimator_output.h"
#include "plate_likelihood_estimator.h"

namespace ilike
{

namespace
{
// a variable's matrix in parameters, or NULL
const arma::mat* find_matrix(const Parameters &parameters, const std::string &variable)
{
  for (auto i = parameters.vector_begin(); i != parameters.vector_end(); ++i)
    if (i->first == variable)
      return i->second.first.get();
  return NULL;
}
}

PlateLikelihoodEstimatorOutput::PlateLikelihoodEstimatorOutput()
  :ExactLikelihoodEstimatorOutput()
{
  this->plate = NULL;
  this->terms_valid = false;
  this->terms_version = 0;
  this->partial = false;
}

PlateLikelihoodEstimatorOutput::PlateLikelihoodEstimatorOutput(PlateLikelihoodEstimator* estimator_in)
  :ExactLikelihoodEstimatorOutput(estimator_in)
{
  this->plate = estimator_in;
  this->terms_valid = false;
  this->terms_version = 0;
  this->partial = false;
}

PlateLikelihoodEstimatorOutput::~PlateLikelihoodEstimatorOutput()
{
}

PlateLikelihoodEstimatorOutput::PlateLikelihoodEstimatorOutput(const PlateLikelihoodEstimatorOutput &another)
  :ExactLikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void PlateLikelihoodEstimatorOutput::operator=(const PlateLikelihoodEstimatorOutput &another)
{
  if (this == &another)
    return;
  ExactLikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* PlateLikelihoodEstimatorOutput::duplicate() const
{
  return new PlateLikelihoodEstimatorOutput(*this);
}

void PlateLikelihoodEstimatorOutput::make_copy(const PlateLikelihoodEstimatorOutput &another)
{
  this->plate = another.plate;
  this->term_values = another.term_values;
  this->terms_valid = another.terms_valid;
  this->terms_version = another.terms_version;
  this->partial = another.partial;
  this->pending = another.pending;
  // (marks are all clear between evaluations)
  this->marked.assign(another.marked.size(), 0);
}

void PlateLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->plate->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->plate_evaluate(parameters);
}

void PlateLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->plate->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->plate_evaluate(parameters);
}

void PlateLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for factors written as model statements.");
}

void PlateLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for factors written as model statements.");
}

arma::mat PlateLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                              const Parameters &x)
{
  Rcpp::stop("Gradients are not available yet for the factor " + this->plate->name +
             " (written as model statements), so gradient-based proposals cannot be used with it.");
}

arma::mat PlateLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                        const Parameters &x)
{
  return this->get_gradient_of_log(variable, x);
}

double PlateLikelihoodEstimatorOutput::plate_evaluate(const Parameters &parameters)
{
  PlateLikelihoodEstimator* p = this->plate;
  size_t n = p->number_of_terms();
  bool use_pending = this->partial && this->terms_valid &&
                     (this->terms_version == p->structure_version) &&
                     (this->term_values.size() == n);
  if (use_pending)
  {
    for (uint32_t t : this->pending)
      this->term_values[t] = p->term(parameters, t);
    if (p->check)
    {
      // check_factor_index: every term must equal its value computed afresh
      for (size_t t = 0; t < n; ++t)
      {
        double fresh = p->term(parameters, t);
        double kept = this->term_values[t];
        if ((fresh != kept) && !(std::isnan(fresh) && std::isnan(kept)) &&
            !(std::abs(fresh - kept) <= 1e-12*(1.0 + std::abs(fresh))))
          Rcpp::stop("check_factor_index: term " + std::to_string(t + 1) + " of the factor " + p->name +
                     " changed although none of the parameters it was listed as reading changed (its list of dependencies is wrong).");
      }
    }
  }
  else
  {
    this->term_values.resize(n);
    for (size_t t = 0; t < n; ++t)
      this->term_values[t] = p->term(parameters, t);
  }
  this->terms_valid = true;
  this->terms_version = p->structure_version;
  this->partial = false;
  this->pending.clear();

  // the sum of all terms in order: the same value as a full evaluation
  double result = 0.0;
  for (size_t t = 0; t < n; ++t)
    result = result + this->term_values[t];
  return result;
}

void PlateLikelihoodEstimatorOutput::carry_over_from(const LikelihoodEstimatorOutput &previous,
                                                     const Parameters &previous_parameters,
                                                     const Parameters &new_parameters)
{
  this->partial = false;
  this->pending.clear();
  this->terms_valid = false;
  const PlateLikelihoodEstimatorOutput* source = dynamic_cast<const PlateLikelihoodEstimatorOutput*>(&previous);
  PlateLikelihoodEstimator* p = this->plate;
  size_t n = p->number_of_terms();
  // terms are carried from an evaluated output of the same plate and
  // structure, between two distinct parameter sets (the comparison below
  // could not see a change made in place to one object)
  if ((source == NULL) || (source == this) || (&previous_parameters == &new_parameters) ||
      (source->plate != p) || !source->terms_valid || source->partial ||
      (source->terms_version != p->structure_version) || (source->term_values.size() != n))
    return;

  if (this->marked.size() != n)
    this->marked.assign(n, 0);
  auto mark = [this](const std::vector<uint32_t> &readers) {
    for (uint32_t t : readers)
      if (!this->marked[t])
      {
        this->marked[t] = 1;
        this->pending.push_back(t);
      }
  };
  bool all = false;
  for (size_t v : p->variables_read)
  {
    const std::string &variable = p->variables[v];
    const arma::mat* before = find_matrix(previous_parameters, variable);
    const arma::mat* after = find_matrix(new_parameters, variable);
    if ((before == NULL) || (after == NULL) ||
        (before->n_rows != after->n_rows) || (before->n_cols != after->n_cols))
    {
      all = true;
      break;
    }
    if ((before == after) || (before->memptr() == after->memptr()))
      continue;
    bool changed = false;
    const auto &elements = p->element_readers[v];
    for (arma::uword k = 0; k < after->n_elem; ++k)
    {
      double a = (*before)(k), b = (*after)(k);
      if ((a == b) || (std::isnan(a) && std::isnan(b)))
        continue;
      changed = true;
      if (!elements.empty())
      {
        uint64_t key = (uint64_t(k % after->n_rows) << 32) | uint64_t(k / after->n_rows);
        auto found = elements.find(key);
        if (found != elements.end())
          mark(found->second);
      }
    }
    if (changed)
      mark(p->whole_readers[v]);
    if (2*this->pending.size() > n)
    {
      // most terms change: evaluate them all
      all = true;
      break;
    }
  }
  for (uint32_t t : this->pending)
    this->marked[t] = 0;
  if (all)
  {
    this->pending.clear();
    return;
  }
  this->term_values = source->term_values;
  this->terms_version = source->terms_version;
  this->terms_valid = true;
  this->partial = true;
}

bool PlateLikelihoodEstimatorOutput::reinitialise(const Parameters &parameters)
{
  if (!ExactLikelihoodEstimatorOutput::reinitialise(parameters))
    return false;
  // the state of a new output: no terms (term_values keeps its memory)
  this->terms_valid = false;
  this->partial = false;
  this->pending.clear();
  return true;
}

bool PlateLikelihoodEstimatorOutput::assign_from(const LikelihoodEstimatorOutput &another)
{
  const PlateLikelihoodEstimatorOutput* source = dynamic_cast<const PlateLikelihoodEstimatorOutput*>(&another);
  if (source == NULL)
    return false;
  LikelihoodEstimatorOutput::make_copy(*source);
  ExactLikelihoodEstimatorOutput::make_copy(*source);
  this->make_copy(*source);
  return true;
}

bool PlateLikelihoodEstimatorOutput::repoint_estimator(LikelihoodEstimator* estimator_in)
{
  PlateLikelihoodEstimator* p = dynamic_cast<PlateLikelihoodEstimator*>(estimator_in);
  if (p == NULL)
    return false;
  this->estimator = p;
  this->plate = p;
  return true;
}

}
