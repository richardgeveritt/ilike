#ifndef DensityLikelihoodEstimatorOutput_H
#define DensityLikelihoodEstimatorOutput_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{
  /**
   * @file density_likelihood_estimator_output.h
   * @brief Defines the DensityLikelihoodEstimator class.
   *
   * Estimates the density likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class DensityLikelihoodEstimator
   * @brief The density likelihood estimator class.
   */


class DensityLikelihoodEstimator;
class DensityEstimator;
class DensityEstimatorOutput;

class DensityLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the densitylikelihoodestimatoroutput operation.
   */
  DensityLikelihoodEstimatorOutput();
  DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
                                   DensityEstimator* density_estimator_in,
                                   DensityEstimator* subsample_density_estimator_in,
                                   bool store_output_in);
  //DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
  //                                 const Parameters &conditioned_on_parameters);
  /**
   * @brief Performs the ~densitylikelihoodestimatoroutput operation.
   */
  virtual ~DensityLikelihoodEstimatorOutput();
  
  /**
   * @brief Performs the densitylikelihoodestimatoroutput operation.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  DensityLikelihoodEstimatorOutput(const DensityLikelihoodEstimatorOutput &another);
  /**
   * @brief Assignment operator for DensityLikelihoodEstimator.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  void operator=(const DensityLikelihoodEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this DensityLikelihoodEstimator object.
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
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void simulate(const Parameters &conditioned_on_parameters);
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
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void subsample_simulate(const Parameters &conditioned_on_parameters);
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
   * @brief Performs the fit operation.
   *
   * @param points The points.
   */
  void fit(const std::vector<Parameters> &points);
  /**
   * @brief Performs the subsample fit operation.
   *
   * @param points The points.
   */
  void subsample_fit(const std::vector<Parameters> &points);
  
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
  DensityLikelihoodEstimator* estimator;
  
  // Stored here.
  /** @brief The density estimator output. */
  DensityEstimatorOutput* density_estimator_output;
  /** @brief The subsample density estimator output. */
  DensityEstimatorOutput* subsample_density_estimator_output;
  // Stored here.
  //DensityEstimator* density_estimator;
  
  //Particles particles;
  /** @brief The points. */
  std::vector<Parameters> points;
  
  /** @brief The time. */
  double time;
  
  /** @brief The log likelihood smcfixed part. */
  double log_likelihood_smcfixed_part;
  /** @brief The subsample log likelihood smcfixed part. */
  double subsample_log_likelihood_smcfixed_part;
  
  /** @brief The store output. */
  bool store_output;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  //std::string results_name;
  
  /** @brief The start time. */
  std::chrono::high_resolution_clock::time_point start_time;
  
  /**
   * @brief Copies the state of another DensityLikelihoodEstimator into this object.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  void make_copy(const DensityLikelihoodEstimatorOutput &another);
  
};
}

#endif
