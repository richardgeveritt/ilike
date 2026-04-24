#ifndef STOCHASTICENSEMBLESHIFTER_H
#define STOCHASTICENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{
  /**
   * @file stochastic_ensemble_shifter.h
   * @brief Defines the StochasticEnsembleShifter class.
   *
   * A stochastic ensemble ensemble shifter. Shifts ensemble members to perform an analysis update in an Ensemble Kalman algorithm.
   *
   * @namespace ilike
   * @class StochasticEnsembleShifter
   * @brief A stochastic ensemble shifter derived from EnsembleShifter.
   */


class StochasticEnsembleShifter : public EnsembleShifter
{
  
public:
  
  /**
   * @brief Default constructor for StochasticEnsembleShifter.
   */
  StochasticEnsembleShifter();
  
  /**
   * @brief Destructor for StochasticEnsembleShifter.
   */
  virtual ~StochasticEnsembleShifter();
  
  /**
   * @brief Copy constructor for StochasticEnsembleShifter.
   *
   * @param another The StochasticEnsembleShifter instance to copy from.
   */
  StochasticEnsembleShifter(const StochasticEnsembleShifter &another);
  
  /**
   * @brief Assignment operator for StochasticEnsembleShifter.
   *
   * @param another The StochasticEnsembleShifter instance to copy from.
   */
  void operator=(const StochasticEnsembleShifter &another);
  /**
   * @brief Creates a deep copy of this StochasticEnsembleShifter object.
   *
   * @return The result.
   */
  EnsembleShifter* duplicate() const;
  
  void setup(Ensemble* ensemble,
             double inverse_incremental_temperature);
  
  void shift(const EnsembleFactorVariables* ensemble_factor_variables,
             arma::colvec &position,
             const std::vector<arma::colvec*> &measurements,
             const std::vector<arma::mat> &kalman_gains,
             const std::vector<arma::colvec> &myys,
             double inverse_incremental_temperature) const;
  
protected:
  
  /**
   * @brief Copies the state of another StochasticEnsembleShifter into this object.
   *
   * @param another The StochasticEnsembleShifter instance to copy from.
   */
  void make_copy(const StochasticEnsembleShifter &another);
  
};
}

#endif
