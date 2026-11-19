#ifndef FACTORVARIABLES_H
#define FACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{

class Particle;
class Index;
class Factors;

class FactorVariables
{

public:

  FactorVariables();
  virtual ~FactorVariables();
  
  FactorVariables(Particle* particle_in);

  FactorVariables(const FactorVariables &another);

  void operator=(const FactorVariables &another);
  virtual FactorVariables* duplicate() const=0;

  // Factor variables for a proposed point new_parameters: the factors in
  // new_outputs_index get new outputs (initialised and simulated at
  // new_parameters, not yet evaluated); every other factor's output is a copy
  // of this object's, keeping its cached value (and, for a random factor, its
  // estimate). Used by MCMC moves, whose other factors either cancel in the
  // acceptance ratio or are refreshed after acceptance. The caller owns the
  // result and must call set_particle on it.
  virtual FactorVariables* copy_with_new_outputs(const Parameters &new_parameters,
                                                 const Index* new_outputs_index) const=0;
  
  virtual void evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  virtual double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  virtual double evaluate_likelihoods(const Index* index) const=0;

  // Sum of the values most recently computed for the factors in index, read
  // from each factor's cache without evaluating anything. Valid only if those
  // factors were last evaluated at the particle's current parameters.
  virtual double get_cached_likelihoods(const Index* index) const=0;

  // Consistency check for MCMC factor indices: for every factor NOT in
  // new_outputs_index whose output is exact (deterministic), re-evaluate a copy
  // at the particle's current parameters and compare with the cached value.
  // Returns the 0-based numbers of the factors whose value changed. Random
  // factors cannot be checked this way and are skipped.
  virtual std::vector<size_t> find_changed_exact_factors(const Index* new_outputs_index) const=0;
  
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  
  virtual double subsample_evaluate_likelihoods(const Index* index) const=0;

  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable) const=0;
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable) const=0;

  // the gradients of several variables at a point other than the particle's
  // (the points of an HMC trajectory, C208): only factors held as a list of
  // outputs (VectorFactorVariables)
  virtual std::vector<arma::mat> direct_get_gradients_of_log_at(const Index* index,
                                                                const std::vector<std::string> &variables,
                                                                const Parameters &at) const
  {
    Rcpp::stop("Gradients at other points are not available for this model's factors (needed by hmc).");
  }

  // the gradients of several variables (C204); by default one at a time
  virtual std::vector<arma::mat> direct_get_gradients_of_log(const Index* index,
                                                             const std::vector<std::string> &variables) const
  {
    std::vector<arma::mat> gradients;
    for (const std::string &variable : variables)
      gradients.push_back(this->direct_get_gradient_of_log(index, variable));
    return gradients;
  }

  virtual const Factors* get_factors() const=0;
  
  void set_particle(Particle* particle_in);
  Particle* get_particle();
  
  virtual void forget_you_were_already_written_to_file()=0;
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index) const=0;

  // Compact records of the factor outputs (see LikelihoodEstimatorOutput::
  // has_compact_value): true if every factor output has one; then
  // compact_values gives one number per factor and write_compact_values
  // writes them as write_to_file would have done. The default (false) keeps
  // full copies of stored iterations.
  // Make this object (the factor variables of a scratch particle) what
  // source->copy_with_new_outputs(new_parameters, new_outputs_index) would
  // return, reusing its output objects. Returns false (and may leave the
  // object partly changed: the caller then makes a new one) if that is not
  // possible. Default: not possible.
  // Make this object what factors->simulate_factor_variables(new_parameters)
  // would return (every output new), reusing its output objects. Returns
  // false if not possible (default).
  virtual bool recycle_all_outputs(const Factors* factors,
                                   const Parameters &new_parameters) { return false; }

  virtual bool recycle_with_new_outputs(const FactorVariables* source,
                                        const Parameters &new_parameters,
                                        const Index* new_outputs_index) { return false; }

  // Adaptive summaries over the factors of an index (dev/plan.md 2.6a, G2):
  // smcadaptive_from_summary returns exactly what
  // evaluate_smcadaptive_part_given_smcfixed_likelihoods(index) would (the
  // sum over the index of each output's value, in index order). Size -1 (the
  // default): not available.
  virtual int adaptive_summary_size(const Index* index) const { return -1; }
  virtual void adaptive_summary(const Index* index, double* values) const {}
  virtual double smcadaptive_from_summary(const Index* index,
                                          const double* values,
                                          const Parameters &parameters) const
  {
    Rcpp::stop("FactorVariables::smcadaptive_from_summary: not available.");
  }
  // the same for n points whose summaries are `stride` numbers apart
  virtual void smcadaptive_from_summaries(const Index* index,
                                          const double* values,
                                          size_t stride,
                                          size_t n,
                                          const Parameters &parameters,
                                          double* results) const
  {
    for (size_t i = 0; i < n; ++i)
      results[i] = this->smcadaptive_from_summary(index, values + i*stride, parameters);
  }

  // Put the outputs in the state of a particle whose factors of index were
  // last evaluated with this adaptive summary (restore_from_summary at
  // parameters) and whose factors of other_index are new (reinitialise);
  // false if not possible (dev/plan.md 2.6a, G3). Default: not possible.
  virtual bool restore_from_summary(const Index* index,
                                    const Index* other_index,
                                    const double* values,
                                    const Parameters &parameters) { return false; }

  // Make these factor variables (and their outputs) belong to `factors`, an
  // object equal to their own (G5); false if not possible (the default).
  virtual bool repoint_to(const Factors* factors) { return false; }

  virtual bool has_compact_values() const { return false; }
  virtual size_t number_of_compact_values() const { return 0; }
  virtual void compact_values(double *values) const {}
  virtual void write_compact_values(const std::string &directory_name,
                                    const std::string &index,
                                    const double *values) const
  {
    Rcpp::stop("FactorVariables::write_compact_values: no compact records.");
  }
  
  virtual void close_ofstreams()=0;

protected:
  
  // not stored here
  Particle* particle;
  
  // not stored here

  void make_copy(const FactorVariables &another);

};
}

#endif
