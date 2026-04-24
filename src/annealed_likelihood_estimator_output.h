#ifndef ANNEALEDLIKELIHOODESTIMATOROUTPUT_H
#define ANNEALEDLIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{
  /**
   * @file annealed_likelihood_estimator_output.h
   * @brief Defines the AnnealedLikelihoodEstimator class.
   *
   * Estimates the annealed likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class AnnealedLikelihoodEstimator
   * @brief The annealed likelihood estimator class.
   */


class AnnealedLikelihoodEstimator;

class AnnealedLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the annealedlikelihoodestimatoroutput operation.
   */
  AnnealedLikelihoodEstimatorOutput();
  AnnealedLikelihoodEstimatorOutput(AnnealedLikelihoodEstimator* estimator_in,
                                    LikelihoodEstimatorOutput* estimator_output_in);
  /**
   * @brief Performs the ~annealedlikelihoodestimatoroutput operation.
   */
  virtual ~AnnealedLikelihoodEstimatorOutput();
  
  /**
   * @brief Performs the annealedlikelihoodestimatoroutput operation.
   *
   * @param another The AnnealedLikelihoodEstimator instance to copy from.
   */
  AnnealedLikelihoodEstimatorOutput(const AnnealedLikelihoodEstimatorOutput &another);
  /**
   * @brief Assignment operator for AnnealedLikelihoodEstimator.
   *
   * @param another The AnnealedLikelihoodEstimator instance to copy from.
   */
  void operator=(const AnnealedLikelihoodEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this AnnealedLikelihoodEstimator object.
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
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  /**
   * @brief Prints the object's state to an output stream.
   *
   * @param os The os.
   */
  void print(std::ostream &os) const;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  /** @brief The estimator. */
  AnnealedLikelihoodEstimator* estimator;
  
  // Stored here.
  /** @brief The estimator output. */
  LikelihoodEstimatorOutput* estimator_output;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  /**
   * @brief Copies the state of another AnnealedLikelihoodEstimator into this object.
   *
   * @param another The AnnealedLikelihoodEstimator instance to copy from.
   */
  void make_copy(const AnnealedLikelihoodEstimatorOutput &another);
  
};
}

#endif
