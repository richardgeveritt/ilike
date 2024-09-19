#ifndef SQUAREROOTENSEMBLESHIFTER_H
#define SQUAREROOTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{
class SquareRootEnsembleShifter : public EnsembleShifter
{
  
public:
  
  SquareRootEnsembleShifter();
  
  virtual ~SquareRootEnsembleShifter();
  
  SquareRootEnsembleShifter(const SquareRootEnsembleShifter &another);
  
  void operator=(const SquareRootEnsembleShifter &another);
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
  
  std::vector<arma::mat> As;
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
  
  void make_copy(const SquareRootEnsembleShifter &another);
  
};
}

#endif
