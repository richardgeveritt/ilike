#ifndef ENSEMBLESHIFTER_H
#define ENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
//#include "ensemble_member.h"

namespace ilike
{
  /**
   * @file ensemble_shifter.h
   * @brief Defines the Ensemble class.
   *
   * Provides ensemble functionality.
   *
   * @namespace ilike
   * @class Ensemble
   * @brief The ensemble class.
   */


class Ensemble;
class EnsembleFactorVariables;

class EnsembleShifter
{
  
public:
  
  /**
   * @brief Performs the ensembleshifter operation.
   */
  EnsembleShifter();
  /**
   * @brief Performs the ~ensembleshifter operation.
   */
  virtual ~EnsembleShifter();
  
  /**
   * @brief Performs the ensembleshifter operation.
   *
   * @param another The Ensemble instance to copy from.
   */
  EnsembleShifter(const EnsembleShifter &another);
  
  /**
   * @brief Assignment operator for Ensemble.
   *
   * @param another The Ensemble instance to copy from.
   */
  void operator=(const EnsembleShifter &another);
  /**
   * @brief Creates a deep copy of this Ensemble object.
   *
   * @return The result.
   */
  virtual EnsembleShifter* duplicate() const=0;
  
  virtual void setup(Ensemble* ensemble,
                     double inverse_incremental_temperature)=0;
  
  virtual void shift(const EnsembleFactorVariables* ensemble_factor_variables,
                     arma::colvec &position,
                     const std::vector<arma::colvec*> &measurements,
                     const std::vector<arma::mat> &kalman_gains,
                     const std::vector<arma::colvec> &myys,
                     double inverse_incremental_temperature) const=0;
  
protected:
  
  /**
   * @brief Copies the state of another Ensemble into this object.
   *
   * @param another The Ensemble instance to copy from.
   */
  void make_copy(const EnsembleShifter &another);
  
};
}

#endif
