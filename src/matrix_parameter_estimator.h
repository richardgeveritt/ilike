#ifndef MATRIXPARAMETERESTIMATOR_H
#define MATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

class MatrixParameterEstimator : public ParameterEstimator
{

public:

  MatrixParameterEstimator();

  virtual ~MatrixParameterEstimator();

  MatrixParameterEstimator(const MatrixParameterEstimator &another);

  void operator=(const MatrixParameterEstimator &another);
  virtual MatrixParameterEstimator* duplicate() const=0;

  virtual void fit(const std::vector<Parameters> &points,
                   arma::colvec normalised_log_weights)=0;
  
  virtual void fit(const std::string &variable,
                   const std::vector<Parameters> &points,
                   arma::colvec normalised_log_weights)=0;
  
  virtual void fit(const std::string &variable,
                   const std::vector<MoveOutput*> &points,
                   arma::colvec normalised_log_weights)=0;
  
  arma::mat estimated;
  
protected:

  void make_copy(const MatrixParameterEstimator &another);

};

#endif
