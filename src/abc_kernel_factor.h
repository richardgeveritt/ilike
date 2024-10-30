#ifndef ABCKERNELFACTOR_H
#define ABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "likelihood_factor.h"
#include "packing_instructions_map.h"

namespace ilike
{

  /**
   * @file abc_kernel_factor.h
   * @brief Defines the ABCKernelFactor class used in ExactLikelihoodEstimator.
   *
   * This file contains the definition of the ABCKernelFactor class, which is a type of LikelihoodFactor
   * used in ExactLikelihoodEstimator. The ABCKernelFactor class has different subtypes corresponding to
   * different types of ABC kernels, each of which is a function of the distance between the observed and
   * simulated data.
   *
   * @namespace ilike
   * @class ABCKernelFactor
   * @brief An ABCKernel type of Factor for use in ExactLikelihoodEstimator.
   *
   * The ABCKernelFactor class is used in ExactLikelihoodEstimator, which contains LikelihoodFactors,
   * DistributionFactors, and IndependentProposalKernels. These are evaluated when a Particle is
   * evaluating its factors in a weight update or Metropolis step, for example.
   */

  class ABCKernelFactor : public LikelihoodFactor
  {

  public:
    /**
     * @brief Default constructor for the ABCKernelFactor class.
     *
     * This constructor initializes an instance of the ABCKernelFactor class.
     */
    ABCKernelFactor();

    /**
     * @brief Constructs an ABCKernelFactor object.
     *
     * @param data_variables_in A vector of strings representing the data variables.
     * @param epsilon_variable_in A string representing the epsilon variable.
     * @param data_in A pointer to a Data object.
     */
    ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                    const std::string &epsilon_variable_in,
                    Data *data_in);

    /**
     * @brief Constructs an ABCKernelFactor object.
     *
     * @param data_variables_in A vector of strings representing the data variables.
     * @param epsilon_variable_in A string representing the epsilon variable.
     * @param scale_variable_in A string representing the scale variable.
     * @param data_in A pointer to a Data object.
     */
    ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                    const std::string &epsilon_variable_in,
                    const std::string &scale_variable_in,
                    Data *data_in);

    /**
     * @brief Virtual destructor for the ABCKernelFactor class.
     *
     * This destructor ensures that derived class destructors are called properly
     * when an object of a derived class is deleted through a pointer to the base class.
     */
    virtual ~ABCKernelFactor();

    /**
     * @brief Copy constructor for the ABCKernelFactor class.
     *
     * This constructor creates a new instance of ABCKernelFactor by copying the
     * data from another existing instance.
     *
     * @param another The ABCKernelFactor instance to copy from.
     */
    ABCKernelFactor(const ABCKernelFactor &another);

    /**
     * @brief Assignment operator for ABCKernelFactor.
     *
     * This operator allows for the assignment of one ABCKernelFactor object to another.
     *
     * @param another The ABCKernelFactor object to be assigned.
     * @return A reference to the assigned ABCKernelFactor object.
     */
    ABCKernelFactor &operator=(const ABCKernelFactor &another);

    /**
     * @brief Creates a deep copy of the current ABCKernelFactor object.
     *
     * This pure virtual function must be implemented by derived classes to
     * provide a mechanism for duplicating the current ABCKernelFactor object.
     *
     * @return A pointer to the duplicated ABCKernelFactor object.
     */
    virtual ABCKernelFactor *abc_kernel_factor_duplicate() const = 0;

    /**
     * @brief Pure virtual function to find the distance based on the given parameters.
     *
     * This function calculates the distance and scale constant based on the input parameters.
     * It must be implemented by any derived class.
     *
     * @param input The parameters used to calculate the distance.
     * @param distance Reference to a double where the calculated distance will be stored.
     * @param scale_constant Reference to a double where the calculated scale constant will be stored.
     */
    virtual void find_distance(const Parameters &input,
                               double &distance,
                               double &scale_constant) const = 0;

    /**
     * @brief Evaluates the kernel function given a distance and a scale constant.
     *
     * @param input The parameters required for the kernel evaluation.
     * @param distance The distance value used in the kernel function.
     * @param scale_constant A constant value used to scale the kernel function.
     * @return The result of the kernel evaluation as a double.
     */
    virtual double evaluate_kernel_given_distance(const Parameters &input,
                                                  double distance,
                                                  double scale_constant) const = 0;

  protected:
    /**
     * @brief Creates a copy of the given ABCKernelFactor object.
     *
     * This function copies the state of the provided ABCKernelFactor
     * instance into the current object.
     *
     * @param another The ABCKernelFactor object to be copied.
     */
    void make_copy(const ABCKernelFactor &another);

    /**
     * @brief A vector that stores a list of data variable names.
     */
    std::vector<std::string> data_variables;

    /**
     * @brief A string variable to store the epsilon value.
     */
    std::string epsilon_variable;

    /**
     * @brief A string representing the scale variable.
     */
    std::string scale_variable;

    /**
     * @brief A column vector from the Armadillo library.
     *
     * A column vector (single column matrix), storing the data.
     */
    arma::colvec data_colvec;

    /**
     * @brief A map that holds packing instructions.
     *
     * This map is used to store and retrieve packing instructions
     * for how to pack data stored in a Data object into a column vector.
     */
    PackingInstructionsMap packing_instructions;
  };
}

#endif
