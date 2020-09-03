#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"



// // From SymbolixAU at https://stackoverflow.com/questions/59055902/find-index-of-all-max-min-values-in-vector-in-rcpp
// IntegerVector which_max(const NumericVector &v)
// {
//   double current_max = v[0];
//   int n = v.size();
//   std::vector< int > res;
//   res.push_back( 0 );
//   int i;
//
//   for( i = 1; i < n; ++i) {
//     double x = v[i];
//     if( x > current_max ) {
//       res.clear();
//       current_max = x;
//       res.push_back( i );
//     } else if ( x == current_max ) {
//       res.push_back( i );
//     }
//   }
//   Rcpp::IntegerVector iv( res.begin(), res.end() );
//   return iv;
// }

double log_sum_exp(const NumericVector &log_weights)
{
  if (sum(is_na(log_weights))>0)
  {
    return R_NegInf;
  }
  unsigned int xmax = which_max(log_weights);
  IntegerVector xmaxv = IntegerVector::create(xmax);
  unsigned int num_weights = log_weights.length();
  IntegerVector all_but_one = setdiff(IntegerVector(seq(0,num_weights-1)),xmaxv);

  NumericVector AllButOnelog_weights = log_weights[all_but_one];

  double result = log1p(sum(exp(AllButOnelog_weights-log_weights[xmax])))+log_weights[xmax];
  if (R_isnancpp(result))
  {
    return R_NegInf;
  }
  else
  {
    return result;
  }
}

NumericVector get_first_element_of_list_as_numeric_vector(const List &a_list)
{
  return a_list[0];
}

// [[Rcpp::export]]
SEXP store_get_first_element_of_list_as_numeric_vector()
{
  return(XPtr<GetDataFromSimulationPtr>(new GetDataFromSimulationPtr(&get_first_element_of_list_as_numeric_vector)));
}

NumericVector identity_statistic(const NumericVector &data)
{
  return data;
}

// [[Rcpp::export]]
SEXP store_identity_statistic()
{
  return(XPtr<SummaryStatisticsPtr>(new SummaryStatisticsPtr(&identity_statistic)));
}

double Lp_uniform_evaluate_log_abc_kernel(const NumericVector &simulated_summary_stats,
                                          const NumericVector &observed_summary_stats,
                                          const double &abc_tolerance,
                                          const double &p)
{
  unsigned int n = simulated_summary_stats.length();
  double Lp_log_distance = (1.0/p)*log_sum_exp(log(abs(simulated_summary_stats-observed_summary_stats)));
  double Lp_ball_log_volume = n*(log(2.0) + lgamma((1.0/p)+1.0) + log(abc_tolerance)) - lgamma((n/p)+1.0);
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume);
  else
    return(R_NegInf);
}

double L1_uniform_evaluate_log_abc_kernel(const NumericVector &simulated_summary_stats,
                                          const NumericVector &observed_summary_stats,
                                          const double &abc_tolerance)
{
  return(Lp_uniform_evaluate_log_abc_kernel(simulated_summary_stats,
                                            observed_summary_stats,
                                            abc_tolerance,
                                            1.0));
}

// [[Rcpp::export]]
SEXP store_L1_uniform_evaluate_log_abc_kernel()
{
  return(XPtr<EvaluateLogABCKernelPtr>(new EvaluateLogABCKernelPtr(&L1_uniform_evaluate_log_abc_kernel)));
}


double L2_uniform_evaluate_log_abc_kernel(const NumericVector &simulated_summary_stats,
                                          const NumericVector &observed_summary_stats,
                                          const double &abc_tolerance)
{
  return(Lp_uniform_evaluate_log_abc_kernel(simulated_summary_stats,
                                            observed_summary_stats,
                                            abc_tolerance,
                                            2.0));
}

// [[Rcpp::export]]
SEXP store_L2_uniform_evaluate_log_abc_kernel()
{
  return(XPtr<EvaluateLogABCKernelPtr>(new EvaluateLogABCKernelPtr(&L2_uniform_evaluate_log_abc_kernel)));
}

double Linf_uniform_evaluate_log_abc_kernel(const NumericVector &simulated_summary_stats,
                                            const NumericVector &observed_summary_stats,
                                            const double &abc_tolerance)
{
  unsigned int n = simulated_summary_stats.length();
  double Lp_log_distance = log(max(abs(simulated_summary_stats-observed_summary_stats)));
  double Lp_ball_log_volume = n*(log(2.0) + log(abc_tolerance));
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume);
  else
    return(R_NegInf);
}

// [[Rcpp::export]]
SEXP store_Linf_uniform_evaluate_log_abc_kernel()
{
  return(XPtr<EvaluateLogABCKernelPtr>(new EvaluateLogABCKernelPtr(&Linf_uniform_evaluate_log_abc_kernel)));
}

double gaussian_evaluate_log_abc_kernel(const NumericVector &simulated_summary_stats,
                                                const NumericVector &observed_summary_stats,
                                                const double &abc_tolerance)
{
  double log_prob = 0.0;
  unsigned int n = simulated_summary_stats.length();
  for (unsigned int i=0; i<n; ++i)
  {
    log_prob = log_prob + R::dnorm(simulated_summary_stats[i]-observed_summary_stats[i],0,abc_tolerance,true);
  }
  return(log_prob);
}

// [[Rcpp::export]]
SEXP store_gaussian_uniform_evaluate_log_abc_kernel()
{
  return(XPtr<EvaluateLogABCKernelPtr>(new EvaluateLogABCKernelPtr(&gaussian_evaluate_log_abc_kernel)));
}

// double EstimateLogLikelihoodUsingEvaluate(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables)
// {
//   return 1;
// }
//
//
// EstimateLogLikelihoodPtr make_estimate_log_likelihood_function_from_evaluate(const EvaluateLogLikelihoodPtr &evaluate_log_likelihood)
// {
//
// }

