//#include <Rcpp.h>
//using namespace Rcpp;

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"

#ifndef ABCLIKELIHOOD_H
#define ABCLIKELIHOOD_H


class ABCLikelihood
{

private:

  double abc_tolerance;

  EvaluateLogABCKernelPtr evaluate_log_abc_kernel;

  SummaryStatisticsPtr summary_statistics;

  //NumericVector summary_statistics_scaling;
  arma::colvec summary_statistics_scaling;

  //NumericVector summary_data;
  arma::colvec summary_data;

protected:

  void MakeCopy(const ABCLikelihood &another);

public:

  ABCLikelihood();

  ABCLikelihood(const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                const SummaryStatisticsPtr &summary_statistics_in,
                const double &abc_tolerance_in,
                const arma::colvec &summary_statistics_scaling_in,
                const List &observed_data_in);

  ABCLikelihood(const ABCLikelihood &another);

  virtual ~ABCLikelihood();

  void operator=(const ABCLikelihood &another);

  double evaluate(const List &simulated_data) const;

  arma::colvec evaluate_multiple(const std::vector<List> &auxiliary_variables) const;

  arma::colvec summary_from_data(const List &simulation) const;

  arma::colvec get_summary_data() const;

  SummaryStatisticsPtr get_summary_statistics() const;

  double get_abc_tolerance() const;

  void set_abc_tolerance(const double &abc_tolerance_in);

  void set_summary_statistics_scaling(const arma::colvec &summary_statistics_scaling_in);

};

#endif
