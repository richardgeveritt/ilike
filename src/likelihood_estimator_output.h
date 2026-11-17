#ifndef LIKELIHOODESTIMATOROUTPUT_H
#define LIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "parameters.h"

namespace ilike
{

class LikelihoodEstimator;
class HMMFactorVariables;
class VectorFactorVariables;
class AnnealedLikelihoodEstimatorOutput;

class LikelihoodEstimatorOutput
{
  
public:
  
  LikelihoodEstimatorOutput();
  
  virtual ~LikelihoodEstimatorOutput();

  LikelihoodEstimatorOutput(const LikelihoodEstimatorOutput &another);
  
  void operator=(const LikelihoodEstimatorOutput &another);
  virtual LikelihoodEstimatorOutput* duplicate() const=0;
  
  // Simulate all auxilliary variables.
  virtual void simulate()=0;
  
  // Simulate all auxilliary variables.
  virtual void simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  double evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  virtual void evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  virtual void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  // Simulate all auxilliary variables.
  virtual void subsample_simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  double subsample_evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  virtual void subsample_evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  virtual LikelihoodEstimator* get_likelihood_estimator() const=0;
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  
  // This will end up a little bit more complicated than we have at the moment, since there is an interaction with
  // calling Simulate.
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Parameters &x)=0;
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Parameters &x)=0;
  
  virtual void forget_you_were_already_written_to_file()=0;
  
  void write(const std::string &directory_name);
  
  friend std::ostream& operator<<(std::ostream& os, const LikelihoodEstimatorOutput &p);
  
  virtual void print(std::ostream &os) const;
  
  double log_likelihood;
  
  double subsample_log_likelihood;
  
  bool write_to_file_flag;
  
  virtual void close_ofstreams()=0;

protected:
  
  friend HMMFactorVariables;
  friend VectorFactorVariables;
  friend AnnealedLikelihoodEstimatorOutput;
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index = "")=0;

  // Compact record of what write_to_file writes, for stored MCMC iterations
  // (dev/plan.md 2.5, "Stage 2 - concrete design", item 3): compact_size()
  // numbers that are enough to write later exactly what write_to_file would
  // write now. 0 (the default) means no compact record: stored iterations
  // then keep a full copy, as before.
  // Reuse of an output object for a new proposal (dev/plan.md, Stage 3 (a)):
  // reinitialise puts this output in exactly the state that
  // estimator->initialise(parameters), then simulate(parameters), then
  // write_to_file_flag = false would give a new one; assign_from copies
  // another output of the same kind (as duplicate() would). Both return false
  // if the output does not support this (the default): a new object is then
  // made, as before.
  virtual bool reinitialise(const Parameters &parameters) { return false; }
  virtual bool assign_from(const LikelihoodEstimatorOutput &another) { return false; }

  virtual size_t compact_size() const { return 0; }

  // Adaptive summaries (dev/plan.md 2.6a, G2). Once the smcfixed part has
  // been evaluated, adaptive_summary() gives adaptive_summary_size() numbers
  // from which smcadaptive_from_summary(values, parameters) computes exactly
  // what evaluate_smcadaptive_part_given_smcfixed(parameters) would set
  // log_likelihood to (same expression), reading from parameters only the
  // SMC sequence variables (annealing power, ABC tolerance). This lets an
  // SMC adapt its next target, or reweight, from a block of numbers instead
  // of the rich particles. smcadaptive_from_summary reads only the summary
  // and the (shared, read-only) estimator, never this output's own state.
  // -1 (the default): not available; the rich path is then used.
  virtual int adaptive_summary_size() const { return -1; }
  virtual void adaptive_summary(double *values) const {}
  virtual double smcadaptive_from_summary(const double *values,
                                          const Parameters &parameters) const
  {
    Rcpp::stop("LikelihoodEstimatorOutput::smcadaptive_from_summary: not available for this output.");
  }

  // smcadaptive_from_summary for n points whose summaries are `stride`
  // numbers apart (results[i] for values + i*stride): the same values, with
  // what does not depend on the point (e.g. the tolerance) looked up once.
  // Default: one call per point.
  virtual void smcadaptive_from_summaries(const double *values,
                                          size_t stride,
                                          size_t n,
                                          const Parameters &parameters,
                                          double *results) const
  {
    for (size_t i = 0; i < n; ++i)
      results[i] = this->smcadaptive_from_summary(values + i*stride, parameters);
  }

  // Put this output in the state its last evaluation left it in (smcfixed
  // part, then the adaptive part at `parameters`), given that evaluation's
  // adaptive summary (dev/plan.md 2.6a, G3). Only for outputs whose
  // estimator supports_restore_from_summary().
  virtual void restore_from_summary(const double *values,
                                    const Parameters &parameters)
  {
    Rcpp::stop("LikelihoodEstimatorOutput::restore_from_summary: not available for this output.");
  }

  // Make this output belong to `estimator` (an estimator equal to its own,
  // e.g. the original of the copy that made it), and the same for every
  // output inside it (dev/plan.md 2.6b, G5). False if not possible (the
  // default).
  virtual bool repoint_estimator(LikelihoodEstimator* estimator_in) { return false; }
  virtual void compact_values(double *values) const {}
  virtual void write_compact_values(const std::string &directory_name,
                                    const std::string &index,
                                    const double *values)
  {
    Rcpp::stop("LikelihoodEstimatorOutput::write_compact_values: this output has no compact record.");
  }
  
  void make_copy(const LikelihoodEstimatorOutput &another);
  
};

}

#endif

