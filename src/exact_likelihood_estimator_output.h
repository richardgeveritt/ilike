#ifndef ExactLikelihoodEstimatorOutput_H
#define ExactLikelihoodEstimatorOutput_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{
  /**
   * @file exact_likelihood_estimator_output.h
   * @brief Defines the ExactLikelihoodEstimator class.
   *
   * Estimates the exact likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class ExactLikelihoodEstimator
   * @brief The exact likelihood estimator class.
   */


class ExactLikelihoodEstimator;

class ExactLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the exactlikelihoodestimatoroutput operation.
   */
  ExactLikelihoodEstimatorOutput();
  /**
   * @brief Performs the exactlikelihoodestimatoroutput operation.
   *
   * @param estimator_in The estimator.
   */
  ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in);
  /**
   * @brief Performs the ~exactlikelihoodestimatoroutput operation.
   */
  virtual ~ExactLikelihoodEstimatorOutput();
  
  /**
   * @brief Performs the exactlikelihoodestimatoroutput operation.
   *
   * @param another The ExactLikelihoodEstimator instance to copy from.
   */
  ExactLikelihoodEstimatorOutput(const ExactLikelihoodEstimatorOutput &another);
  /**
   * @brief Assignment operator for ExactLikelihoodEstimator.
   *
   * @param another The ExactLikelihoodEstimator instance to copy from.
   */
  void operator=(const ExactLikelihoodEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this ExactLikelihoodEstimator object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* duplicate() const;
  
  /**
   * @brief Simulates the required variables.
   */
  void simulate();
  
  /**
   * @brief Simulates the required variables.
   *
   * @param parameters The parameters.
   */
  void simulate(const Parameters &parameters);
  //double evaluate(const Parameters &parameters);
  
  /**
   * @brief Evaluates the smcfixed part.
   *
   * @param parameters The parameters.
   */
  void evaluate_smcfixed_part(const Parameters &parameters);
  /**
   * @brief Evaluates the smcadaptive part given smcfixed.
   *
   * @param parameters The parameters.
   */
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  //void evaluate_smcfixed_part(const Parameters &parameters);
  //void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  /**
   * @brief Performs the subsample simulate operation.
   */
  void subsample_simulate();
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param parameters The parameters.
   */
  void subsample_simulate(const Parameters &parameters);
  
  /**
   * @brief Performs the subsample evaluate smcfixed part operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed operation.
   *
   * @param parameters The parameters.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  /**
   * @brief Returns the likelihood estimator.
   *
   * @return The result.
   */
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /**
   * @brief Prints the object's state to an output stream.
   *
   * @param os The os.
   */
  void print(std::ostream &os) const;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  /** @brief The estimator. */
  ExactLikelihoodEstimator* estimator;
  
  /** @brief The log likelihood smcfixed part. */
  double log_likelihood_smcfixed_part;
  /** @brief The subsample log likelihood smcfixed part. */
  double subsample_log_likelihood_smcfixed_part;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  /**
   * @brief Copies the state of another ExactLikelihoodEstimator into this object.
   *
   * @param another The ExactLikelihoodEstimator instance to copy from.
   */
  void make_copy(const ExactLikelihoodEstimatorOutput &another);
  
};
}

#endif
