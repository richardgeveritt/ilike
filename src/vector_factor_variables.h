#ifndef VECTORFACTORVARIABLES_H
#define VECTORFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <memory>
#include <string>
#include "factor_variables.h"

namespace ilike
{

class VectorFactors;

class LikelihoodEstimatorOutput;

class VectorFactorVariables : public FactorVariables
{
  
public:
  
  VectorFactorVariables();
  
  VectorFactorVariables(const VectorFactors* vector_factors,
                        std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in);

  VectorFactorVariables(VectorFactors* vector_factors,
                        std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs_in,
                        Particle* particle_in);
  
  virtual ~VectorFactorVariables();
  
  VectorFactorVariables(const VectorFactorVariables &another);
  
  void operator=(const VectorFactorVariables &another);
  FactorVariables* duplicate() const;
  FactorVariables* copy_with_new_outputs(const Parameters &new_parameters,
                                         const Index* new_outputs_index) const;
  
  void evaluate_smcfixed_part_of_likelihoods(const Index* index);
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  
  double evaluate_likelihoods(const Index* index) const;
  double get_cached_likelihoods(const Index* index) const;
  std::vector<size_t> find_changed_exact_factors(const Index* new_outputs_index) const;
  
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);

  double subsample_evaluate_likelihoods(const Index* index) const;
  
  arma::mat direct_get_gradient_of_log(const std::string &variable) const;
  
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable) const;
  
  std::vector<arma::mat> direct_get_gradients_of_log(const Index* index,
                                                     const std::vector<std::string> &variables) const;
  arma::mat direct_get_gradient_of_log(const Index* index,
                                       const std::string &variable) const;
  
  arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                 const std::string &variable) const;
  
  const Factors* get_factors() const;
  
  void forget_you_were_already_written_to_file();
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;

  bool recycle_with_new_outputs(const FactorVariables* source,
                                const Parameters &new_parameters,
                                const Index* new_outputs_index);

  bool recycle_all_outputs(const Factors* factors,
                           const Parameters &new_parameters);

  int adaptive_summary_size(const Index* index) const;
  void adaptive_summary(const Index* index, double* values) const;
  double smcadaptive_from_summary(const Index* index,
                                  const double* values,
                                  const Parameters &parameters) const;
  void smcadaptive_from_summaries(const Index* index,
                                  const double* values,
                                  size_t stride,
                                  size_t n,
                                  const Parameters &parameters,
                                  double* results) const;

  bool restore_from_summary(const Index* index,
                            const Index* other_index,
                            const double* values,
                            const Parameters &parameters);

  bool repoint_to(const Factors* factors);

  bool has_compact_values() const;
  size_t number_of_compact_values() const;
  void compact_values(double *values) const;
  void write_compact_values(const std::string &directory_name,
                            const std::string &index,
                            const double *values) const;
  
  void close_ofstreams();
  
  // In-place moves (dev/plan.md I4c-2): factor i's output; a new output of
  // factor i for parameters (made and simulated as a proposal's would be),
  // put in place of the current one, which is returned.
  LikelihoodEstimatorOutput* get_output(size_t i) const;
  std::unique_ptr<LikelihoodEstimatorOutput> replace_with_new_output(size_t i, const Parameters &parameters);
  void put_back_output(size_t i, std::unique_ptr<LikelihoodEstimatorOutput> output);

protected:
  
  // not stored here
  const VectorFactors* vector_factors;
  
  void make_copy(const VectorFactorVariables &another);
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<LikelihoodEstimatorOutput>> likelihood_estimator_outputs;

  // The factor directory names (directory + "/factorK") for the directory
  // they were last built for: write_compact_values is called once per stored
  // iteration with the same directory, and building the names each time was
  // most of the cost of writing a chain's factor records. Not copied. Behind
  // a pointer (NULL until first used) so that the many particles that never
  // write compact records (IS, SMC) carry 8 bytes, not two containers.
  struct FactorDirectories
  {
    std::string directory;
    std::vector<std::string> names;
  };
  mutable std::unique_ptr<FactorDirectories> factor_directories;

};
}

#endif
