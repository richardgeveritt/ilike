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

/*
double log_sum_exp(const NumericVector &log_weights)
{
  if (sum(is_na(log_weights))>0)
  {
    return R_NegInf;
  }
  //unsigned int xmax = which_max(log_weights);
  //IntegerVector xmaxv = IntegerVector::create(xmax);
  //unsigned int num_weights = log_weights.length();
  //IntegerVector all_but_one = setdiff(IntegerVector(seq(0,num_weights-1)),xmaxv);
  //NumericVector AllButOnelog_weights = log_weights[all_but_one];
  //double result = log1p(sum(exp(AllButOnelog_weights-log_weights[xmax])))+log_weights[xmax];
  double xmax = max(log_weights);
  double result = xmax + log(sum(exp(log_weights-xmax)));

  if (R_isnancpp(result))
  {
    return R_NegInf;
  }
  else
  {
    return result;
  }

}
*/

double log_sum_exp(const std::vector<double> &log_weights)
{
  double xmax = log_weights[0];
  for (std::vector<double>::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  {
    if ((*i)>xmax)
      xmax = *i;
  }
  double the_sum = 0;
  for (std::vector<double>::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  {
    the_sum = the_sum + exp((*i)-xmax);
  }
  return(log(the_sum));
}

double log_sum_exp(const arma::colvec &log_weights)
{
  double xmax = log_weights[0];
  for (arma::colvec::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  {
    if ((*i)>xmax)
      xmax = *i;
  }
  double the_sum = 0;
  for (arma::colvec::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  {
    the_sum = the_sum + exp((*i)-xmax);
  }
  return(log(the_sum));
}

std::vector<double> operator*(const std::vector<double> &first,
                              const std::vector<double> &second)
{
  std::vector<double> output;
  std::size_t n = first.size();
  output.reserve(n);
  for (unsigned int i =0; i<n; ++i)
  {
    output.push_back(first[i]*second[i]);
  }
  return output;
}

std::vector<double> operator-(const std::vector<double> &first,
                              const std::vector<double> &second)
{
  std::vector<double> output;
  std::size_t n = first.size();
  output.reserve(n);
  for (unsigned int i =0; i<n; ++i)
  {
    output.push_back(first[i]-second[i]);
  }
  return output;
}

std::vector<double> abs(const std::vector<double> &first)
{
  std::vector<double> output;
  std::size_t n = first.size();
  output.reserve(n);
  for (std::vector<double>::const_iterator i = first.begin(); i!=first.end(); ++i)
  {
    output.push_back(abs(*i));
  }
  return(output);
}

std::vector<double> log(const std::vector<double> &first)
{
  std::vector<double> output;
  std::size_t n = first.size();
  output.reserve(n);
  for (std::vector<double>::const_iterator i = first.begin(); i!=first.end(); ++i)
  {
    output.push_back(log(*i));
  }
  return(output);
}

std::vector<double> operator*(const double &c,
                              const std::vector<double> &first)
{
  std::vector<double> output;
  std::size_t n = first.size();
  output.reserve(n);
  for (std::vector<double>::const_iterator i = first.begin(); i!=first.end(); ++i)
  {
    output.push_back(c*(*i));
  }
  return(output);
}

double cess(const arma::colvec &normalised_log_weights,
            const arma::colvec &log_incremental_weights)
{
  if ( sum(log_incremental_weights==R_NegInf)==normalised_log_weights.size() )
  {
    return(0);
  }
  else if ( sum(log_incremental_weights==R_PosInf)==normalised_log_weights.size() )
  {
    return(log_incremental_weights.size());
  }
  else
  {
    return(exp(log(log_incremental_weights.size()) + 2.0*(log_sum_exp(normalised_log_weights + log_incremental_weights)) - log_sum_exp(normalised_log_weights + 2.0*log_incremental_weights)));
  }
}

NumericMatrix get_first_element_of_list_as_numeric_matrix(const List &a_list)
{
  return a_list[0];
}

// [[Rcpp::export]]
SEXP store_get_first_element_of_list_as_numeric_matrix()
{
  return(XPtr<GetDataFromSimulationPtr>(new GetDataFromSimulationPtr(&get_first_element_of_list_as_numeric_matrix)));
}

arma::colvec make_vector_statistic(const NumericMatrix &data)
{
  arma::colvec output(data.length());
  unsigned int counter = 0;
  for (unsigned int j=0; j<data.ncol(); ++j)
  {
    for (unsigned int i=0; i<data.nrow(); ++i)
    {
      output[counter] = data(i,j);
      counter = counter + 1;
    }
  }
  return output;
}

// [[Rcpp::export]]
SEXP store_make_vector_statistic()
{
  return(XPtr<SummaryStatisticsPtr>(new SummaryStatisticsPtr(&make_vector_statistic)));
}

double Lp_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
                                          const double &abc_tolerance,
                                          const double &p)
{
  unsigned int n = simulated_summary_stats.size();
  double Lp_log_distance = (1.0/p)*log_sum_exp(p*log(abs(simulated_summary_stats-observed_summary_stats)));
  double Lp_ball_log_volume = double(n)*(log(2.0) + lgamma((1.0/p)+1.0) + log(abc_tolerance)) - lgamma((double(n)/p)+1.0);
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume);
  else
    return(R_NegInf);
}

double L1_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
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


double L2_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
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

double Linf_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                            const arma::colvec &observed_summary_stats,
                                            const double &abc_tolerance)
{
  unsigned int n = simulated_summary_stats.size();
  double Lp_log_distance = log(max(abs(simulated_summary_stats-observed_summary_stats)));
  double Lp_ball_log_volume = double(n)*(log(2.0) + log(abc_tolerance));
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

double gaussian_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                        const arma::colvec &observed_summary_stats,
                                        const double &abc_tolerance)
{
  double log_prob = 0.0;
  unsigned int n = simulated_summary_stats.size();
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

