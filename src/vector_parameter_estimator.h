#ifndef VECTORPARAMETERESTIMATOR_H
#define VECTORPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

namespace ilike
{
class VectorParameterEstimator : public ParameterEstimator
{
  
public:
  
  VectorParameterEstimator();
  
  virtual ~VectorParameterEstimator();
  
  VectorParameterEstimator(const VectorParameterEstimator &another);
  
  void operator=(const VectorParameterEstimator &another);
  virtual VectorParameterEstimator* duplicate() const=0;
  
  virtual void fit(const arma::mat &points,
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
  
  arma::colvec estimated;
  
protected:
  
  void make_copy(const VectorParameterEstimator &another);
  
};
}

#endif
