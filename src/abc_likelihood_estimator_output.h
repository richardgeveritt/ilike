#ifndef ABCLIKELIHOODESTIMATOROUTPUT_H
#define ABCLIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{

  /**
   * @file abc_likelihood_estimator_output.h
   * @brief Header file for the ABCLikelihoodEstimatorOutput class.
   *
   * This file contains the declaration of the ABCLikelihoodEstimatorOutput class,
   * which is a subclass of LikelihoodEstimatorOutput. It provides methods for
   * simulating and evaluating likelihoods, both in full and subsample contexts.
   */

  /**
   * @class ABCLikelihoodEstimatorOutput
   * @brief A class for handling output from an ABC likelihood estimator.
   *
   * This class provides methods to simulate and evaluate likelihoods,
   * manage gradient calculations, and handle file operations related to
   * the output of an ABC likelihood estimator.
   */

  class ABCLikelihoodEstimator;

  class ABCLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
  {

  public:
    /**
     * @brief Default constructor.
     */
    ABCLikelihoodEstimatorOutput();

    /**
     * @brief Constructor with an estimator.
     * @param estimator_in Pointer to an ABCLikelihoodEstimator.
     */
    ABCLikelihoodEstimatorOutput(ABCLikelihoodEstimator *estimator_in);

    /**
     * @brief Destructor.
     */
    virtual ~ABCLikelihoodEstimatorOutput();

    /**
     * @brief Copy constructor.
     * @param another Another instance of ABCLikelihoodEstimatorOutput to copy from.
     */
    ABCLikelihoodEstimatorOutput(const ABCLikelihoodEstimatorOutput &another);

    /**
     * @brief Assignment operator.
     * @param another Another instance of ABCLikelihoodEstimatorOutput to assign from.
     */
    void operator=(const ABCLikelihoodEstimatorOutput &another);

    /**
     * @brief Creates a deep copy of the current object.
     * @return A pointer to the duplicated LikelihoodEstimatorOutput.
     */
    LikelihoodEstimatorOutput *duplicate() const;

    /**
     * @brief Simulates the variables used in the estimator.
     */
    void simulate();

    /**
     * @brief Simulates the variables used in the estimator given the parameters.
     * @param parameters The parameters to use for simulation.
     */
    void simulate(const Parameters &parameters);

    /**
     * @brief Evaluates the part of the kernel that does not change at each iteration of the SMC.
     *
     * @param parameters A constant reference to a Parameters object containing the necessary data for evaluation.
     */
    void evaluate_smcfixed_part(const Parameters &parameters);

    /**
     * @brief Evaluates the kernel that does change at each iteration of the SMC.
     *
     * This function performs the evaluation of the kernel part
     * using the provided parameters. It is used in the context of SMC algorithms where
     * certain parameters are fixed and others are adapted during the process.
     *
     * @param parameters The fixed parameters used for the evaluation.
     */
    void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

    /**
     * @brief Simulates the variables used in the estimator, in the case where a subsample of data is used, given parameters.
     * @param parameters The parameters to use for subsample simulation.
     */
    void subsample_simulate(const Parameters &parameters);

    /**
     * @brief Evaluates the part of the kernel that does not change at each iteration of the SMC, in the case where the data is subsampled.
     *
     * @param parameters A constant reference to a Parameters object containing the necessary data for evaluation.
     */
    void subsample_evaluate_smcfixed_part(const Parameters &parameters);

    /**
     * @brief Evaluates the part of the kernel that does not change at each iteration of the SMC, in the case where the data is subsampled, given parameters.
     *
     * @param parameters A constant reference to a Parameters object containing the necessary data for evaluation.
     */
    void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

    /**
     * @brief Retrieves the associated LikelihoodEstimator object.
     *
     * This function returns a pointer to the LikelihoodEstimator object
     * associated with this instance. The returned pointer is constant,
     * indicating that the LikelihoodEstimator cannot be modified through
     * this pointer.
     *
     * @return LikelihoodEstimator* A constant pointer to the associated
     * LikelihoodEstimator object.
     */
    LikelihoodEstimator *get_likelihood_estimator() const;

    /**
     * @brief Computes the gradient of the log-likelihood with respect to a given variable.
     *
     * @param variable The name of the variable with respect to which the gradient is computed.
     * @param x The parameters used in the computation of the gradient.
     * @return arma::mat The gradient of the log-likelihood as a matrix.
     */
    arma::mat get_gradient_of_log(const std::string &variable,
                                  const Parameters &x);

    /**
     * @brief Computes the gradient of the log-likelihood with respect to a given variable, when the data is subsampled.
     *
     * @param variable The name of the variable with respect to which the gradient is computed.
     * @param x The parameters used in the computation of the gradient.
     * @return arma::mat The gradient of the log-likelihood as a matrix.
     */
    arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                            const Parameters &x);

    /**
     * @brief Resets the state indicating that the object has been written to a file.
     *
     * This function clears any internal flags or states that track whether the
     * object has already been written to a file, allowing it to be written again
     * if necessary.
     */
    void forget_you_were_already_written_to_file();

    /**
     * @brief Closes all open output file streams.
     *
     * This function ensures that any output file streams that were opened
     * are properly closed, releasing any associated resources.
     */
    void close_ofstreams();

    /**
     * @brief Prints the contents of the object to the given output stream.
     *
     * This function outputs the internal state or relevant information of the object
     * to the provided output stream. It is useful for debugging or logging purposes.
     *
     * @param os The output stream where the object's information will be printed.
     */
    void print(std::ostream &os) const;

  protected:
    // Stored in sampler.
    /**
     * @brief Pointer to an ABCLikelihoodEstimator instance.
     *
     * This member variable holds a pointer to an ABCLikelihoodEstimator object,
     * which is responsible for performing likelihood estimation using Approximate
     * Bayesian Computation (ABC) methods.
     */
    ABCLikelihoodEstimator *estimator;

    /**
     * @brief Stores the distance.
     */
    double distance;

    /**
     * @brief A constant used to scale statistics within the likelihood estimator.
     *
     * This variable is used to adjust the scale of the different statistics in the likelihood estimator.
     */
    double scale_constant;

    /**
     * @brief Writes the output to a file in the specified directory.
     *
     * @param directory_name The name of the directory where the file will be written.
     * @param index An optional index to append to the file name. Default is an empty string.
     */
    void write_to_file(const std::string &directory_name,
                       const std::string &index = "");

    /**
     * @brief Creates a copy of the given ABCLikelihoodEstimatorOutput object.
     *
     * This function duplicates the state of the provided ABCLikelihoodEstimatorOutput
     * instance into the current instance.
     *
     * @param another The ABCLikelihoodEstimatorOutput object to be copied.
     */
    void make_copy(const ABCLikelihoodEstimatorOutput &another);
  };
}

#endif
