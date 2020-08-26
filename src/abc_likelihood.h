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

  NumericVector tolerances;

  EvaluateLogABCKernelPtr evaluate_log_abc_kernel;

  SummaryStatisticPtr summary_statistic;

  NumericVector summary_statistic_scaling;

  NumericVector summary_data;

protected:

  void MakeCopy(const ABCLikelihood &another);

public:

  ABCLikelihood();

  ABCLikelihood(const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                const SummaryStatisticPtr &summary_statistic_in,
                const NumericVector &tolerances_in,
                const NumericVector &summary_statistic_scaling_in,
                const NumericVector &data_in);

  ABCLikelihood(const ABCLikelihood &another);

  virtual ~ABCLikelihood();

  void operator=(const ABCLikelihood &another);

  double evaluate(const NumericVector &simulated_data) const;

  NumericVector get_summary_data() const;

  SummaryStatisticPtr get_summary_statistic() const;

  void set_tolerances(const NumericVector &tolerances_in);

  void set_summary_statistic_scaling(const NumericVector &summary_statistic_scaling_in);

};

#endif
