#include <Rcpp.h>
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

  NumericVector summary_statistics_scaling;

  NumericVector summary_data;

protected:

  void MakeCopy(const ABCLikelihood &another);

public:

  ABCLikelihood();

  ABCLikelihood(const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                const SummaryStatisticsPtr &summary_statistics_in,
                const double &abc_tolerance_in,
                const NumericVector &summary_statistics_scaling_in,
                const NumericVector &data_in);

  ABCLikelihood(const ABCLikelihood &another);

  virtual ~ABCLikelihood();

  void operator=(const ABCLikelihood &another);

  double evaluate(const NumericVector &simulated_data) const;

  NumericVector get_summary_data() const;

  SummaryStatisticsPtr get_summary_statistics() const;

  void set_abc_tolerance(const double &abc_tolerance_in);

  void set_summary_statistics_scaling(const NumericVector &summary_statistics_scaling_in);

};

#endif
