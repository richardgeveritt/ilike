#ifndef LIKELIHOODESTIMATOROUTPUT_H
#define LIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "parameters.h"

namespace ilike
{
  /**
   * @file likelihood_estimator_output.h
   * @brief Defines the LikelihoodEstimator class.
   *
   * Estimates the likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class LikelihoodEstimator
   * @brief The likelihood estimator class.
   */


class LikelihoodEstimator;
class HMMFactorVariables;
class VectorFactorVariables;
class AnnealedLikelihoodEstimatorOutput;
//class Data;

class LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the likelihoodestimatoroutput operation.
   */
  LikelihoodEstimatorOutput();
  
  /**
   * @brief Performs the ~likelihoodestimatoroutput operation.
   */
  virtual ~LikelihoodEstimatorOutput();
  
  //virtual void simulate(const Parameters &parameters);
  
  /**
   * @brief Performs the likelihoodestimatoroutput operation.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  LikelihoodEstimatorOutput(const LikelihoodEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for LikelihoodEstimator.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  void operator=(const LikelihoodEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this LikelihoodEstimator object.
   *
   * @return The result.
   */
  virtual LikelihoodEstimatorOutput* duplicate() const=0;
  
  // Simulate all auxilliary variables.
  /**
   * @brief Simulates the required variables.
   */
  virtual void simulate()=0;
  
  // Simulate all auxilliary variables.
  /**
   * @brief Simulates the required variables.
   *
   * @param parameters The parameters.
   */
  virtual void simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  /**
   * @brief Evaluates.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  /**
   * @brief Evaluates the smcfixed part.
   *
   * @param parameters The parameters.
   */
  virtual void evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  /**
   * @brief Evaluates the smcadaptive part given smcfixed.
   *
   * @param parameters The parameters.
   */
  virtual void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  // Simulate all auxilliary variables.
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param parameters The parameters.
   */
  virtual void subsample_simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  /**
   * @brief Performs the subsample evaluate operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double subsample_evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  /**
   * @brief Performs the subsample evaluate smcfixed part operation.
   *
   * @param parameters The parameters.
   */
  virtual void subsample_evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed operation.
   *
   * @param parameters The parameters.
   */
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  /**
   * @brief Returns the likelihood estimator.
   *
   * @return The result.
   */
  virtual LikelihoodEstimator* get_likelihood_estimator() const=0;
  
  /**
   * @brief Performs the change data operation.
   */
  void change_data();
  /**
   * @brief Performs the change data operation.
   *
   * @param new_data The new data.
   */
  void change_data(std::shared_ptr<Data> new_data);
  
  // This will end up a little bit more complicated than we have at the moment, since there is an interaction with
  // calling Simulate.
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Parameters &x)=0;
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Parameters &x)=0;
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  virtual void forget_you_were_already_written_to_file()=0;
  
  /**
   * @brief Writes results to a file.
   *
   * @param directory_name The directory name.
   */
  void write(const std::string &directory_name);
  
  /**
   * @brief Performs the operator<< operation.
   *
   * @param os The os.
   * @param p The p.
   *
   * @return The result.
   */
  friend std::ostream& operator<<(std::ostream& os, const LikelihoodEstimatorOutput &p);
  
  /**
   * @brief Prints the object's state to an output stream.
   *
   * @param os The os.
   */
  virtual void print(std::ostream &os) const;
  
  double log_likelihood;
  
  double subsample_log_likelihood;
  
  bool write_to_file_flag;
  
  /**
   * @brief Closes any open file streams.
   */
  virtual void close_ofstreams()=0;
  
  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  //
  // virtual List simulate_auxiliary_variables(const List &inputs) const=0;
  //
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
  
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

  // Put this output in the state its last evaluation left it in (smcfixed
  // part, then the adaptive part at `parameters`), given that evaluation's
  // adaptive summary (dev/plan.md 2.6a, G3). Only for outputs whose
  // estimator supports_restore_from_summary().
  virtual void restore_from_summary(const double *values,
                                    const Parameters &parameters)
  {
    Rcpp::stop("LikelihoodEstimatorOutput::restore_from_summary: not available for this output.");
  }
  virtual void compact_values(double *values) const {}
  virtual void write_compact_values(const std::string &directory_name,
                                    const std::string &index,
                                    const double *values)
  {
    Rcpp::stop("LikelihoodEstimatorOutput::write_compact_values: this output has no compact record.");
  }
  
  /**
   * @brief Copies the state of another LikelihoodEstimator into this object.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  void make_copy(const LikelihoodEstimatorOutput &another);
  
};

}

#endif

