#include "plate_likelihood_estimator_output.h"
#include "plate_likelihood_estimator.h"

namespace ilike
{

namespace
{
// a variable's matrix in parameters, or NULL
const arma::mat* find_matrix(const Parameters &parameters, const std::string &variable)
{
  return parameters.find_vector(variable);
}
}

PlateLikelihoodEstimatorOutput::PlateLikelihoodEstimatorOutput()
  :ExactLikelihoodEstimatorOutput()
{
  this->plate = NULL;
  this->terms_valid = false;
  this->terms_version = 0;
  this->partial = false;
  this->leaves = 1;
  this->undo_log_likelihood = 0.0;
  this->undo_smcfixed_part = 0.0;
  this->undo_pending = false;
}

PlateLikelihoodEstimatorOutput::PlateLikelihoodEstimatorOutput(PlateLikelihoodEstimator* estimator_in)
  :ExactLikelihoodEstimatorOutput(estimator_in)
{
  this->plate = estimator_in;
  this->terms_valid = false;
  this->terms_version = 0;
  this->partial = false;
  this->leaves = 1;
  this->undo_log_likelihood = 0.0;
  this->undo_smcfixed_part = 0.0;
  this->undo_pending = false;
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
  this->tree = another.tree;
  this->leaves = another.leaves;
  this->undo = another.undo;
  this->undo_log_likelihood = another.undo_log_likelihood;
  this->undo_smcfixed_part = another.undo_smcfixed_part;
  this->undo_pending = another.undo_pending;
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
  if (this->plate->gradient_function == NULL)
    Rcpp::stop("Gradients are not available for the factor " + this->plate->name +
               " (written as model statements), so gradient-based proposals cannot be used with it.");
  return this->plate->gradient_function(variable, x, this->plate->plate_data());
}

// several variables in one pass over the terms (C204)
std::vector<arma::mat> PlateLikelihoodEstimatorOutput::get_gradients_of_log(const std::vector<std::string> &variables,
                                                                            const Parameters &x)
{
  if (this->plate->gradients_function == NULL)
    return LikelihoodEstimatorOutput::get_gradients_of_log(variables, x);
  std::vector<arma::mat> gradients;
  this->plate->gradients_function(variables, gradients, x, this->plate->plate_data());
  return gradients;
}

bool PlateLikelihoodEstimatorOutput::evaluate_with_gradients(const Parameters &parameters,
                                                             const std::vector<std::string> &variables,
                                                             double &value,
                                                             std::vector<arma::mat> &gradients)
{
  // (check_factor_index checks the terms, which this does not keep)
  if ((this->plate->value_gradients_function == NULL) || this->plate->check)
    return false;
  gradients.clear();
  value = this->plate->value_gradients_function(variables, gradients, parameters, this->plate->plate_data());
  this->undo.clear();
  this->undo_pending = false;
  this->terms_valid = false;
  this->partial = false;
  this->pending.clear();
  this->set_value(value);
  return true;
}

arma::mat PlateLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                        const Parameters &x)
{
  return this->get_gradient_of_log(variable, x);
}

void PlateLikelihoodEstimatorOutput::size_tree(size_t n)
{
  size_t m = 1;
  while (m < n)
    m = 2*m;
  this->leaves = m;
  this->tree.assign(2*m, 0.0);
}

void PlateLikelihoodEstimatorOutput::build_tree()
{
  for (size_t i = this->leaves - 1; i >= 1; --i)
    this->tree[i] = this->tree[2*i] + this->tree[2*i + 1];
}

void PlateLikelihoodEstimatorOutput::update_tree(uint32_t t)
{
  for (size_t i = (this->leaves + t)/2; i >= 1; i = i/2)
    this->tree[i] = this->tree[2*i] + this->tree[2*i + 1];
}

double PlateLikelihoodEstimatorOutput::tree_sum() const
{
  // (with one leaf, node 1 is that leaf)
  return this->tree[1];
}

// the tree holds n terms
static bool tree_fits(const std::vector<double> &tree, size_t leaves, size_t n)
{
  size_t m = 1;
  while (m < n)
    m = 2*m;
  return (leaves == m) && (tree.size() == 2*m);
}

void PlateLikelihoodEstimatorOutput::set_value(double value)
{
  if (this->plate->smcfixed_flag)
    this->log_likelihood_smcfixed_part = value;
  this->log_likelihood = value;
}

double PlateLikelihoodEstimatorOutput::plate_evaluate(const Parameters &parameters)
{
  PlateLikelihoodEstimator* p = this->plate;
  size_t n = p->number_of_terms();
  this->undo.clear();
  this->undo_pending = false;
  bool use_pending = this->partial && this->terms_valid &&
                     (this->terms_version == p->structure_version) &&
                     tree_fits(this->tree, this->leaves, n);
  if (use_pending)
  {
    for (uint32_t t : this->pending)
      this->tree[this->leaves + t] = p->term(parameters, t);
    if (4*this->pending.size() > n)
      this->build_tree();
    else
      for (uint32_t t : this->pending)
        this->update_tree(t);
    if (p->check)
    {
      // check_factor_index: every term must equal its value computed afresh
      for (size_t t = 0; t < n; ++t)
      {
        double fresh = p->term(parameters, t);
        double kept = this->tree[this->leaves + t];
        if ((fresh != kept) && !(std::isnan(fresh) && std::isnan(kept)) &&
            !(std::abs(fresh - kept) <= 1e-12*(1.0 + std::abs(fresh))))
          Rcpp::stop("check_factor_index: term " + std::to_string(t + 1) + " of the factor " + p->name +
                     " changed although none of the parameters it was listed as reading changed (its list of dependencies is wrong).");
      }
    }
  }
  else
  {
    this->size_tree(n);
    for (size_t t = 0; t < n; ++t)
      this->tree[this->leaves + t] = p->term(parameters, t);
    this->build_tree();
  }
  this->terms_valid = true;
  this->terms_version = p->structure_version;
  this->partial = false;
  this->pending.clear();
  return this->tree_sum();
}

bool PlateLikelihoodEstimatorOutput::change_in_place(const Parameters &parameters,
                                                     const std::string &variable,
                                                     size_t row,
                                                     size_t column)
{
  PlateLikelihoodEstimator* p = this->plate;
  size_t n = p->number_of_terms();
  if (!this->terms_valid || this->partial || (this->terms_version != p->structure_version) ||
      !tree_fits(this->tree, this->leaves, n))
    return false;
  this->undo.clear();
  this->undo_log_likelihood = this->log_likelihood;
  this->undo_smcfixed_part = this->log_likelihood_smcfixed_part;
  this->undo_pending = true;
  size_t v = p->variables.size();
  for (size_t k = 0; k < p->variables.size(); ++k)
    if (p->variables[k] == variable)
      v = k;
  if (v == p->variables.size())
    return true; // no term reads the variable: nothing changes

  if (this->marked.size() != n)
    this->marked.assign(n, 0);
  auto recompute = [&](const std::vector<uint32_t> &readers) {
    for (uint32_t t : readers)
      if (!this->marked[t])
      {
        this->marked[t] = 1;
        this->undo.push_back(std::make_pair(t, this->tree[this->leaves + t]));
        this->tree[this->leaves + t] = p->term(parameters, t);
      }
  };
  const auto &elements = p->element_readers[v];
  auto found = elements.find((uint64_t(row) << 32) | uint64_t(column));
  if (found != elements.end())
    recompute(found->second);
  recompute(p->whole_readers[v]);
  for (const auto &changed : this->undo)
  {
    this->marked[changed.first] = 0;
    this->update_tree(changed.first);
  }
  if (p->check)
  {
    for (size_t t = 0; t < n; ++t)
    {
      double fresh = p->term(parameters, t);
      double kept = this->tree[this->leaves + t];
      if ((fresh != kept) && !(std::isnan(fresh) && std::isnan(kept)) &&
          !(std::abs(fresh - kept) <= 1e-12*(1.0 + std::abs(fresh))))
        Rcpp::stop("check_factor_index: term " + std::to_string(t + 1) + " of the factor " + p->name +
                   " changed although none of the parameters it was listed as reading changed (its list of dependencies is wrong).");
    }
  }
  this->set_value(this->tree_sum());
  return true;
}

void PlateLikelihoodEstimatorOutput::undo_in_place()
{
  if (!this->undo_pending)
    return;
  for (auto i = this->undo.rbegin(); i != this->undo.rend(); ++i)
    this->tree[this->leaves + i->first] = i->second;
  for (const auto &changed : this->undo)
    this->update_tree(changed.first);
  this->log_likelihood = this->undo_log_likelihood;
  this->log_likelihood_smcfixed_part = this->undo_smcfixed_part;
  this->undo.clear();
  this->undo_pending = false;
}

void PlateLikelihoodEstimatorOutput::prepare_in_place(const Parameters &parameters)
{
  PlateLikelihoodEstimator* p = this->plate;
  if (this->terms_valid && !this->partial && (this->terms_version == p->structure_version) &&
      tree_fits(this->tree, this->leaves, p->number_of_terms()))
    return;
  // (the terms at the current parameters, whose value is already known: the
  // tree is made valid, the value kept)
  this->plate_evaluate(parameters);
}

void PlateLikelihoodEstimatorOutput::commit_in_place()
{
  this->undo.clear();
  this->undo_pending = false;
}

void PlateLikelihoodEstimatorOutput::carry_over_from(const LikelihoodEstimatorOutput &previous,
                                                     const Parameters &previous_parameters,
                                                     const Parameters &new_parameters,
                                                     const std::vector<std::string>* changed)
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
      (source->terms_version != p->structure_version) || !tree_fits(source->tree, source->leaves, n))
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
  // compare the values of variable v before and after, and mark the terms
  // reading its changed elements; true: give up (evaluate every term)
  auto examine = [&](size_t v) -> bool {
    const std::string &variable = p->variables[v];
    const arma::mat* before = find_matrix(previous_parameters, variable);
    const arma::mat* after = find_matrix(new_parameters, variable);
    if ((before == NULL) || (after == NULL) ||
        (before->n_rows != after->n_rows) || (before->n_cols != after->n_cols))
      return true;
    if ((before == after) || (before->memptr() == after->memptr()))
      return false;
    bool differs = false;
    const auto &elements = p->element_readers[v];
    for (arma::uword k = 0; k < after->n_elem; ++k)
    {
      double a = (*before)(k), b = (*after)(k);
      if ((a == b) || (std::isnan(a) && std::isnan(b)))
        continue;
      differs = true;
      if (!elements.empty())
      {
        uint64_t key = (uint64_t(k % after->n_rows) << 32) | uint64_t(k / after->n_rows);
        auto found = elements.find(key);
        if (found != elements.end())
          mark(found->second);
      }
    }
    if (differs)
      mark(p->whole_readers[v]);
    // most terms change: evaluate them all
    return (2*this->pending.size() > n);
  };
  bool all = false;
  if (changed != NULL)
  {
    // only the variables the move changed (C205): the others are the same
    // matrices, so their terms cannot have changed
    for (const std::string &name : *changed)
    {
      auto found = p->read_index.find(name);
      if ((found != p->read_index.end()) && examine(found->second))
      {
        all = true;
        break;
      }
    }
  }
  else
  {
    for (size_t v : p->variables_read)
      if (examine(v))
      {
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
  this->tree = source->tree;
  this->leaves = source->leaves;
  this->terms_version = source->terms_version;
  this->terms_valid = true;
  this->partial = true;
}

bool PlateLikelihoodEstimatorOutput::reinitialise(const Parameters &parameters)
{
  if (!ExactLikelihoodEstimatorOutput::reinitialise(parameters))
    return false;
  // the state of a new output: no terms (the tree keeps its memory)
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
