#ifndef ANNEALEDLIKELIHOODESTIMATOROUTPUT_H
#define ANNEALEDLIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{

class AnnealedLikelihoodEstimator;

class AnnealedLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  AnnealedLikelihoodEstimatorOutput();
  AnnealedLikelihoodEstimatorOutput(AnnealedLikelihoodEstimator* estimator_in,
                                    LikelihoodEstimatorOutput* estimator_output_in);
  virtual ~AnnealedLikelihoodEstimatorOutput();
  
  AnnealedLikelihoodEstimatorOutput(const AnnealedLikelihoodEstimatorOutput &another);
  void operator=(const AnnealedLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  
  void simulate();
  void simulate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  std::vector<arma::mat> get_gradients_of_log(const std::vector<std::string> &variables,
                                              const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  void close_ofstreams();

  // passed on to the annealed output (dev/plan.md I4b)
  void carry_over_from(const LikelihoodEstimatorOutput &previous,
                       const Parameters &previous_parameters,
                       const Parameters &new_parameters,
                       const std::vector<std::string>* changed = NULL);
  
  // in-place moves on one element, passed on to the annealed output; the
  // value is the power times its new value (dev/plan.md I4c-2)
  bool change_in_place(const Parameters &parameters, const std::string &variable,
                       size_t row, size_t column);
  void undo_in_place();
  void commit_in_place();
  void prepare_in_place(const Parameters &parameters);

  void forget_you_were_already_written_to_file();
  
  void print(std::ostream &os) const;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  AnnealedLikelihoodEstimator* estimator;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<LikelihoodEstimatorOutput> estimator_output;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");

  size_t compact_size() const;
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;
  void restore_from_summary(const double *values,
                            const Parameters &parameters);
  bool repoint_estimator(LikelihoodEstimator* estimator_in);
  
  void compact_values(double *values) const;
  
  void write_compact_values(const std::string &directory_name,
  
                            const std::string &index,
  
                            const double *values);
  
  void make_copy(const AnnealedLikelihoodEstimatorOutput &another);

  // the power at the parameters (as in evaluate_smcadaptive_part_given_smcfixed)
  double power(const Parameters &parameters) const;

  // the value before an in-place change
  double undo_log_likelihood = 0.0;
  
};
}

#endif
