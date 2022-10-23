#ifndef STOCHASTICENSEMBLESHIFTER_H
#define STOCHASTICENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

class StochasticEnsembleShifter : public EnsembleShifter
{

public:

  StochasticEnsembleShifter();

  virtual ~StochasticEnsembleShifter();

  StochasticEnsembleShifter(const StochasticEnsembleShifter &another);

  void operator=(const StochasticEnsembleShifter &another);
  EnsembleShifter* duplicate() const;
  
  void setup(Ensemble* ensemble);
  
  void shift(const EnsembleFactorVariables* ensemble_factor_variables,
             arma::colvec &position,
             const std::vector<arma::mat> &Cxys,
             const std::vector<arma::mat> &Cyys,
             double inverse_incremental_temperature) const;
  
protected:

  void make_copy(const StochasticEnsembleShifter &another);

};

#endif
