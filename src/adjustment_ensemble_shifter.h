#ifndef ADJUSTMENTENSEMBLESHIFTER_H
#define ADJUSTMENTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{

class AdjustmentEnsembleShifter : public EnsembleShifter
{
  
public:
  
  AdjustmentEnsembleShifter();
  
  virtual ~AdjustmentEnsembleShifter();
  
  AdjustmentEnsembleShifter(const AdjustmentEnsembleShifter &another);
  
  void operator=(const AdjustmentEnsembleShifter &another);
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
  
  arma::mat Zf;
  //arma::mat Ginv;
  //arma::mat Ftranspose;
  arma::colvec mean_position;
  //std::vector<arma::mat> Vs;
  std::vector<arma::mat> As;
  
  arma::mat P;
  arma::mat Vtranspose;
  arma::mat Dhathalf;
  std::vector<arma::mat> Yhats;
  
  void make_copy(const AdjustmentEnsembleShifter &another);
  
};
}

#endif
