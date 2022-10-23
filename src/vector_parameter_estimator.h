#ifndef VECTORPARAMETERESTIMATOR_H
#define VECTORPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

class VectorParameterEstimator : public ParameterEstimator
{

public:

  VectorParameterEstimator();

  virtual ~VectorParameterEstimator();

  VectorParameterEstimator(const VectorParameterEstimator &another);

  void operator=(const VectorParameterEstimator &another);
  virtual VectorParameterEstimator* duplicate() const=0;

  virtual void fit(const std::vector<Parameters> &points,
                   arma::colvec normalised_log_weights)=0;
  
  virtual void fit(const std::string &variable,
                   const std::vector<Parameters> &points,
                   arma::colvec normalised_log_weights)=0;
  
  virtual void fit(const std::string &variable,
                   const std::vector<MoveOutput*> &points,
                   arma::colvec normalised_log_weights)=0;
  
  arma::colvec estimated;
  
protected:

  void make_copy(const VectorParameterEstimator &another);

};

#endif
