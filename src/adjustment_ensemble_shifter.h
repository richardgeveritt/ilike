#ifndef ADJUSTMENTENSEMBLESHIFTER_H
#define ADJUSTMENTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{
  /**
   * @file adjustment_ensemble_shifter.h
   * @brief Defines the AdjustmentEnsembleShifter class.
   *
   * A adjustment ensemble ensemble shifter. Shifts ensemble members to perform an analysis update in an Ensemble Kalman algorithm.
   *
   * @namespace ilike
   * @class AdjustmentEnsembleShifter
   * @brief An adjustment ensemble shifter derived from EnsembleShifter.
   */



class AdjustmentEnsembleShifter : public EnsembleShifter
{
  
public:
  
  /**
   * @brief Default constructor for AdjustmentEnsembleShifter.
   */
  AdjustmentEnsembleShifter();
  
  /**
   * @brief Destructor for AdjustmentEnsembleShifter.
   */
  virtual ~AdjustmentEnsembleShifter();
  
  /**
   * @brief Copy constructor for AdjustmentEnsembleShifter.
   *
   * @param another The AdjustmentEnsembleShifter instance to copy from.
   */
  AdjustmentEnsembleShifter(const AdjustmentEnsembleShifter &another);
  
  /**
   * @brief Assignment operator for AdjustmentEnsembleShifter.
   *
   * @param another The AdjustmentEnsembleShifter instance to copy from.
   */
  void operator=(const AdjustmentEnsembleShifter &another);
  /**
   * @brief Creates a deep copy of this AdjustmentEnsembleShifter object.
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
  
  /** @brief The zf. */
  arma::mat Zf;
  //arma::mat Ginv;
  //arma::mat Ftranspose;
  /** @brief The mean position. */
  arma::colvec mean_position;
  //std::vector<arma::mat> Vs;
  /** @brief The as. */
  std::vector<arma::mat> As;
  
  /** @brief The p. */
  arma::mat P;
  /** @brief The vtranspose. */
  arma::mat Vtranspose;
  /** @brief The dhathalf. */
  arma::mat Dhathalf;
  /** @brief The yhats. */
  std::vector<arma::mat> Yhats;
  
  /**
   * @brief Copies the state of another AdjustmentEnsembleShifter into this object.
   *
   * @param another The AdjustmentEnsembleShifter instance to copy from.
   */
  void make_copy(const AdjustmentEnsembleShifter &another);
  
};
}

#endif
