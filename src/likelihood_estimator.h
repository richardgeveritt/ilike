#include <Rcpp.h>
using namespace Rcpp;

#include <vector>

#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

class LikelihoodEstimator
{

protected:

  NumericVector data;

public:

  LikelihoodEstimator(const NumericVector &data_in);

  virtual ~LikelihoodEstimator();

  virtual double estimate_log_likelihood(const NumericVector &inputs,
                                         const List &auxiliary_variables) const=0;

  virtual List simulate_auxiliary_variables(const NumericVector &inputs) const=0;

  virtual void setup_likelihood_estimator(const NumericMatrix &all_points,
                                          const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
