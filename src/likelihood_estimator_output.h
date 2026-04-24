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
  
  /**
   * @brief Copies the state of another LikelihoodEstimator into this object.
   *
   * @param another The LikelihoodEstimator instance to copy from.
   */
  void make_copy(const LikelihoodEstimatorOutput &another);
  
};

}

#endif

