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

