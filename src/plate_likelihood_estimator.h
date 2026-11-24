#ifndef PLATELIKELIHOODESTIMATOR_H
#define PLATELIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include "exact_likelihood_estimator.h"

namespace ilike
{

// A factor made of many terms, each reading a few elements of the parameters
// (dev/plan.md 3.6 Layer 1, I4b). Written by the ilike language for model
// statements: the structure function lists the terms (one row each: the
// statement, then its loop indices) and what each term reads (rows: term,
// variable, row, column, with row = column = -1 for the whole variable,
// 0-based, variables numbered in `variables`); the term function gives one
// term's log density. The factor's value is the sum of its terms, in order.
// Each output keeps its terms' values, so that a move recomputes only the
// terms whose inputs changed (PlateLikelihoodEstimatorOutput).
typedef Data (*PlateStructurePtr)(const Data &data);
// (dev/plan.md 3.21 (c5)) a factor that is one ABC kernel: a point's two
// summary numbers, and the kernel at the current tolerance from them
typedef void (*AbcSummaryPtr)(const Parameters &parameters, const Data &data, double *values);
typedef double (*AbcFromSummaryPtr)(const double *values, const Parameters &parameters, const Data &data);
typedef double (*PlateTermPtr)(const Parameters &parameters,
                               const Data &data,
                               const arma::mat &terms,
                               size_t term);
// the factor's gradient with respect to one variable (generated from its
// terms, dev/plan.md I4e)
typedef arma::mat (*PlateGradientPtr)(const std::string &variable,
                                      const Parameters &parameters,
                                      const Data &data);

// the gradients of several variables in one pass over the terms (C204)
typedef void (*PlateGradientsPtr)(const std::vector<std::string> &variables,
                                  std::vector<arma::mat> &gradients,
                                  const Parameters &parameters,
                                  const Data &data);

// the factor's value and the gradients of several variables in one pass (C222)
typedef double (*PlateValueGradientsPtr)(const std::vector<std::string> &variables,
                                         std::vector<arma::mat> &gradients,
                                         const Parameters &parameters,
                                         const Data &data);

class PlateLikelihoodEstimatorOutput;

class PlateLikelihoodEstimator : public ExactLikelihoodEstimator
{
public:

  PlateLikelihoodEstimator();
  PlateLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           PlateStructurePtr structure_in,
                           PlateTermPtr term_in,
                           const std::vector<std::string> &variables_in,
                           const std::string &name_in,
                           bool check_in,
                           bool smcfixed_flag_in);
  virtual ~PlateLikelihoodEstimator();

  PlateLikelihoodEstimator(const PlateLikelihoodEstimator &another);
  void operator=(const PlateLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);

  size_t number_of_terms() const { return this->terms.n_rows; }

  // data tempering (dev/plan.md 3.20, S2): the terms are put in tempering
  // order (a random order fixed by the seed, or the order given) when the
  // structure is built; at a value phi of the sequence variable, with n terms
  // and n phi = k + r, the target has the first k terms and term k + 1 to the
  // power r (anneal), or the first ceil(n phi) terms
  void set_data_tempering(const std::string &variable_in, bool anneal_in, bool random_order_in, size_t seed_in);
  bool is_tempered() const { return this->tempered; }
  void tempering_point(const Parameters &parameters, size_t &k, double &r) const;
  void set_gradient(PlateGradientPtr gradient_in) { this->gradient_function = gradient_in; }
  // (dev/plan.md 3.21 (c4)) the factor's largest value (generated from its
  // terms' maxima), for early rejection; without it, not known (+Inf). (It
  // used to be ExactLikelihoodEstimator's, which with no factors of its own
  // gave 0, as if known.)
  void set_maximum(EvaluateLogLikelihoodPtr maximum_in) { this->maximum_function = maximum_in; }
  // (dev/plan.md 3.21 (c5)) an ABC kernel's summaries: with them the factor
  // supports flat generations even when it changes along the SMC sequence
  // (an outer tolerance), its points weighted from two numbers each
  void set_kernel_summary(AbcSummaryPtr summary_in, AbcFromSummaryPtr from_summary_in)
  {
    this->summary_function = summary_in;
    this->from_summary_function = from_summary_in;
  }
  bool has_kernel_summary() const { return (this->summary_function != NULL) && !this->tempered; }
  void kernel_summary(const Parameters &parameters, double *values) const { this->summary_function(parameters, *this->data, values); }
  double kernel_from_summary(const double *values, const Parameters &parameters) const
  {
    return this->from_summary_function(values, parameters, *this->data);
  }
  bool supports_flat_generations() const { return this->smcfixed_flag || this->has_kernel_summary(); }
  double log_maximum(const Parameters &parameters) const
  {
    if (this->maximum_function == NULL)
      return arma::datum::inf;
    return this->maximum_function(parameters, *this->data);
  }
  void set_gradients(PlateGradientsPtr gradients_in) { this->gradients_function = gradients_in; }
  void set_value_gradients(PlateValueGradientsPtr value_gradients_in) { this->value_gradients_function = value_gradients_in; }

protected:
  AbcSummaryPtr summary_function = NULL;
  AbcFromSummaryPtr from_summary_function = NULL;
  EvaluateLogLikelihoodPtr maximum_function = NULL;

  friend PlateLikelihoodEstimatorOutput;

  void specific_change_data(Data* new_data);

  // (re)reads the structure from the data; increments structure_version
  void build_structure();
  const Data &plate_data() const;
  double term(const Parameters &parameters, size_t t) const;

  PlateStructurePtr structure;
  PlateTermPtr term_function;
  PlateGradientPtr gradient_function = NULL;
  PlateGradientsPtr gradients_function = NULL;
  PlateValueGradientsPtr value_gradients_function = NULL;
  std::vector<std::string> variables;
  std::string name;
  // check_factor_index: after a partial evaluation, recompute every term and
  // stop if one not recomputed changed (a wrong dependency list)
  bool check;

  // data tempering (set_data_tempering)
  bool tempered = false;
  std::string tempering_variable;
  bool tempering_anneal = true;
  bool tempering_random = true;
  size_t tempering_seed = 0;

  Data* current_plate_data;
  bool structure_built;
  // outputs carry term values only within one version of the structure
  size_t structure_version;
  arma::mat terms;
  // the terms reading each variable as a whole, and each element of it
  // (key: row << 32 | column)
  std::vector<std::vector<uint32_t>> whole_readers;
  std::vector<std::unordered_map<uint64_t, std::vector<uint32_t>>> element_readers;
  // the variables some term reads, and their numbers by name (C205)
  std::vector<size_t> variables_read;
  std::unordered_map<std::string, size_t> read_index;

  void make_copy(const PlateLikelihoodEstimator &another);
};

}

#endif
