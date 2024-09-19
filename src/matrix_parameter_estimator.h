#ifndef MATRIXPARAMETERESTIMATOR_H
#define MATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

namespace ilike
{
class MatrixParameterEstimator : public ParameterEstimator
{
  
public:
  
  MatrixParameterEstimator();
  
  virtual ~MatrixParameterEstimator();
  
  MatrixParameterEstimator(const MatrixParameterEstimator &another);
  
  void operator=(const MatrixParameterEstimator &another);
  virtual MatrixParameterEstimator* duplicate() const=0;
  
  virtual void fit(const arma::mat &matrix_points,
                   const arma::colvec &normalised_log_weights)=0;
  
  void fit(const std::vector<std::string> &variables,
           const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<MoveOutput*> &points,
           const arma::colvec &normalised_log_weights);
  
  arma::mat estimated;
  
protected:
  
  void make_copy(const MatrixParameterEstimator &another);
  
};
}

#endif
