#ifndef SQUAREROOTENSEMBLESHIFTER_H
#define SQUAREROOTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{
  /**
   * @file square_root_ensemble_shifter.h
   * @brief Defines the SquareRootEnsembleShifter class.
   *
   * A square root ensemble ensemble shifter. Shifts ensemble members to perform an analysis update in an Ensemble Kalman algorithm.
   *
   * @namespace ilike
   * @class SquareRootEnsembleShifter
   * @brief A square root ensemble shifter derived from EnsembleShifter.
   */


class SquareRootEnsembleShifter : public EnsembleShifter
{
  
public:
  
  /**
   * @brief Default constructor for SquareRootEnsembleShifter.
   */
  SquareRootEnsembleShifter();
  
  /**
   * @brief Destructor for SquareRootEnsembleShifter.
   */
  virtual ~SquareRootEnsembleShifter();
  
  /**
   * @brief Copy constructor for SquareRootEnsembleShifter.
   *
   * @param another The SquareRootEnsembleShifter instance to copy from.
   */
  SquareRootEnsembleShifter(const SquareRootEnsembleShifter &another);
  
  /**
   * @brief Assignment operator for SquareRootEnsembleShifter.
   *
   * @param another The SquareRootEnsembleShifter instance to copy from.
   */
  void operator=(const SquareRootEnsembleShifter &another);
  /**
   * @brief Creates a deep copy of this SquareRootEnsembleShifter object.
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
  
  /** @brief The as. */
  std::vector<arma::mat> As;
  /** @brief The mean position. */
  arma::colvec mean_position;
  
  //arma::mat Zf;
  //arma::mat Ginv;
  //arma::mat Ftranspose;
  //arma::colvec mean_position;
  //std::vector<arma::mat> Vs;
  //arma::mat Sigma;
  //std::vector<arma::mat> HSigmaHts;
  
  //arma::mat P;
  //arma::mat Vtranspose;
  //arma::mat Dhathalf;
  //std::vector<arma::mat> Yhats;
  
  /**
   * @brief Copies the state of another SquareRootEnsembleShifter into this object.
   *
   * @param another The SquareRootEnsembleShifter instance to copy from.
   */
  void make_copy(const SquareRootEnsembleShifter &another);
  
};
}

#endif
