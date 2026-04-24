#ifndef ZEROFINDINGDOUBLEPARAMETERESTIMATOR_H
#define ZEROFINDINGDOUBLEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "double_recursive_parameter_estimator.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file zero_finding_double_recursive_parameter_estimator.h
   * @brief Defines the ZeroFindingDoubleRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for zero finding double parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class ZeroFindingDoubleRecursiveParameterEstimator
   * @brief A zero finding double recursive parameter estimator derived from DoubleRecursiveParameterEstimator.
   */


class ZeroFindingDoubleRecursiveParameterEstimator : public DoubleRecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for ZeroFindingDoubleRecursiveParameterEstimator.
   */
  ZeroFindingDoubleRecursiveParameterEstimator();
  ZeroFindingDoubleRecursiveParameterEstimator(double initial_value,
                                               double target_score_in);
  
  /**
   * @brief Destructor for ZeroFindingDoubleRecursiveParameterEstimator.
   */
  virtual ~ZeroFindingDoubleRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for ZeroFindingDoubleRecursiveParameterEstimator.
   *
   * @param another The ZeroFindingDoubleRecursiveParameterEstimator instance to copy from.
   */
  ZeroFindingDoubleRecursiveParameterEstimator(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for ZeroFindingDoubleRecursiveParameterEstimator.
   *
   * @param another The ZeroFindingDoubleRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this ZeroFindingDoubleRecursiveParameterEstimator object.
   *
   * @return The result.
   */
  RecursiveParameterEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a double pointer.
   *
   * @return The result.
   */
  DoubleRecursiveParameterEstimator* double_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  /** @brief The gain. */
  GainPtr gain;
  
  /** @brief The target score. */
  double target_score;
  
  /**
   * @brief Copies the state of another ZeroFindingDoubleRecursiveParameterEstimator into this object.
   *
   * @param another The ZeroFindingDoubleRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const ZeroFindingDoubleRecursiveParameterEstimator &another);
  
};
}

#endif
