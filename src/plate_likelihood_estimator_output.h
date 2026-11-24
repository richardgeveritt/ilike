#ifndef PLATELIKELIHOODESTIMATOROUTPUT_H
#define PLATELIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <cstdint>
#include <vector>
#include "exact_likelihood_estimator_output.h"

namespace ilike
{

class PlateLikelihoodEstimator;

// The output of a plate factor (dev/plan.md 3.6 Layer 1, I4b): the value of
// each term at this particle's parameters. When a move makes the output of a
// proposed particle, carry_over_from copies the current particle's terms and
// marks the terms that read an element whose value differs between the two
// parameter sets; the next evaluation recomputes only those. Every other way
// of reaching a state (reinitialise, restoring from a summary or a compact
// record) leaves no valid terms, so the next evaluation recomputes them all.
// The factor's value is always the sum of all the terms in order, so a
// partial evaluation gives exactly the value of a full one.
class PlateLikelihoodEstimatorOutput : public ExactLikelihoodEstimatorOutput
{
public:

  PlateLikelihoodEstimatorOutput();
  PlateLikelihoodEstimatorOutput(PlateLikelihoodEstimator* estimator_in);
  virtual ~PlateLikelihoodEstimatorOutput();

  PlateLikelihoodEstimatorOutput(const PlateLikelihoodEstimatorOutput &another);
  void operator=(const PlateLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  // (dev/plan.md 3.21 (c5)) an ABC kernel factor that changes along the
  // sequence: its two summary numbers, kept when the smcfixed part is evaluated
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values, const Parameters &parameters) const;
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  std::vector<arma::mat> get_gradients_of_log(const std::vector<std::string> &variables,
                                              const Parameters &x);
  // (C222) the value by one pass over the terms with the gradients: the
  // terms themselves are not kept, so the output is left with no valid terms
  // (the next evaluation, or a scan's prepare_in_place, recomputes them)
  bool evaluate_with_gradients(const Parameters &parameters,
                               const std::vector<std::string> &variables,
                               double &value,
                               std::vector<arma::mat> &gradients);

  bool carries_over() const { return true; }
  void carry_over_from(const LikelihoodEstimatorOutput &previous,
                       const Parameters &previous_parameters,
                       const Parameters &new_parameters,
                       const std::vector<std::string>* changed = NULL);

  // In-place moves on one element (I4c-2): after parameters[variable](row,
  // column) was changed in place, recompute the terms reading that element,
  // keeping their old values; false if this output has no valid terms (the
  // caller then evaluates a new output instead). undo_in_place restores the
  // terms and the value; commit_in_place keeps the change.
  bool change_in_place(const Parameters &parameters,
                       const std::string &variable,
                       size_t row,
                       size_t column);
  void undo_in_place();
  void commit_in_place();
  void prepare_in_place(const Parameters &parameters);

protected:
  double kernel_values[2] = {0.0, 0.0};

  bool reinitialise(const Parameters &parameters);
  bool assign_from(const LikelihoodEstimatorOutput &another);
  bool repoint_estimator(LikelihoodEstimator* estimator_in);

  double plate_evaluate(const Parameters &parameters);

  // data tempering (dev/plan.md 3.20, S2): only the terms the target
  // includes are evaluated; the first `active` terms are valid (the others are
  // not, and are never read). The value is the sum of the first k terms plus
  // r times term k + 1.
  size_t active = 0;
  double tempered_evaluate(const Parameters &parameters);
  double tempered_value(const Parameters &parameters) const;
  double prefix_sum(size_t k) const;

  // Stored in ModelAndAlgorithm.
  PlateLikelihoodEstimator* plate;

  // The terms' values as the leaves of a summation tree: leaves - 1 inner
  // nodes, then the leaves (term t at tree[leaves + t], zeros after the last
  // term); node i is the sum of nodes 2i and 2i + 1, and node 1 is the
  // factor's value. Updating k terms costs O(k log n), and a full evaluation
  // builds the same tree, so both give exactly the same value. The terms
  // are those at this output's parameters (except those listed in pending,
  // when partial), for the plate's structure_version.
  std::vector<double> tree;
  size_t leaves;
  void size_tree(size_t n);
  void build_tree();
  void update_tree(uint32_t t);
  double tree_sum() const;
  void set_value(double value);
  // in-place changes not yet committed: (term, old value), and the old value
  std::vector<std::pair<uint32_t, double>> undo;
  double undo_log_likelihood;
  double undo_smcfixed_part;
  bool undo_pending;
  bool terms_valid;
  size_t terms_version;
  bool partial;
  std::vector<uint32_t> pending;
  std::vector<char> marked;

  void make_copy(const PlateLikelihoodEstimatorOutput &another);
};

}

#endif
