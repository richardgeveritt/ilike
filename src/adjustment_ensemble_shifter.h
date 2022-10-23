#ifndef ADJUSTMENTENSEMBLESHIFTER_H
#define ADJUSTMENTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

class AdjustmentEnsembleShifter : public EnsembleShifter
{

public:

  AdjustmentEnsembleShifter();
  
  virtual ~AdjustmentEnsembleShifter();

  AdjustmentEnsembleShifter(const AdjustmentEnsembleShifter &another);

  void operator=(const AdjustmentEnsembleShifter &another);
  EnsembleShifter* duplicate() const;
  
  void setup(Ensemble* ensemble);
  
  void shift(const EnsembleFactorVariables* ensemble_factor_variables,
             arma::colvec &position,
             const std::vector<arma::mat> &Cxys,
             const std::vector<arma::mat> &Cyys,
             double inverse_incremental_temperature) const;

protected:
  
  arma::mat Zf;
  arma::mat Ginv;
  arma::mat Ftranspose;
  arma::colvec mean_position;
  std::vector<arma::mat> Vs;

  void make_copy(const AdjustmentEnsembleShifter &another);

};

#endif
