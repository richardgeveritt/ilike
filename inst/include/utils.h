//#include <Rcpp.h>
//using namespace Rcpp;

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <functional>
//#include "ilike_header.h"
#include "parameters.h"
//#include <RcppCommon.h>



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

//double log_sum_exp(const NumericVector &log_weights);

//double log_sum_exp(const std::vector<double> &log_weights);

inline double log_sum_exp(const arma::colvec &log_weights)
{
  //double xmax = log_weights[0];
  //for (arma::colvec::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  //{
  //  if ((*i)>xmax)
  //    xmax = *i;
  //}
  
  size_t which_max = log_weights.index_max();
  double xmax = log_weights(which_max);
  if (xmax==-arma::datum::inf)
    return -arma::datum::inf;
  //double the_sum = 0;
  //for (arma::colvec::const_iterator i = log_weights.begin(); i!=log_weights.end(); ++i)
  //{
  //  the_sum = the_sum + exp((*i)-xmax);
  //}
  arma::colvec log_weights_copy = log_weights;
  log_weights_copy.shed_row(which_max);
  
  double the_sum = arma::sum(arma::exp(log_weights_copy-xmax));
  return(xmax+log1p(the_sum));
}

inline std::vector<size_t> stratified_resample(const arma::colvec &log_weights,
                                               const arma::colvec &resampling_variables)
{
  size_t n = log_weights.size();
  std::vector<size_t> indices(n);
  
  arma::colvec norm_log_weights = log_weights - log_sum_exp(log_weights);
  arma::colvec W = exp(norm_log_weights);
  arma::colvec cw = cumsum(W / sum(W));
  //std::cout << cw.t() << std::endl;
  arma::colvec resampling_variables_over_n = resampling_variables/double(n);
  
  arma::colvec seq(n);
  for (size_t i=0; i<n; ++i)
  {
    seq[i] = i;
  }
  
  arma::colvec v = resampling_variables_over_n + seq/double(n);
  
  //std::cout << v.t() << std::endl;
  
  //std::cout << resampling_variables_over_n.t() << std::endl;

  size_t j = 0;
  for (size_t i=0; i<n; ++i)
  {
    while(cw[j] < v[i])
    {
      j = j+1;
    }
    indices[i] = j;
  }

  return(indices);
}

inline arma::rowvec mean_wt(const arma::mat &x,
                            const arma::colvec &wt)
{
  if (wt.n_rows!=x.n_rows)
    Rcpp::stop("meanwt - matrix and weight sizes don't match.");
  arma::mat wt_mat = repmat(wt,1,x.n_cols);
  return arma::sum(wt_mat % x,0);
}

inline arma::mat cov_wt(const arma::mat &x,
                        const arma::colvec &wt)
{
  arma::colvec mean = mean_wt(x,wt);
  arma::mat x_minus_mean = x-repmat(mean,x.n_rows,1);
  arma::mat wt_mat = repmat(wt,1,x.n_cols);
  return ((wt_mat % x_minus_mean).t()*x_minus_mean)/(1 - sum(wt%wt));
}

inline arma::mat vector_of_parameters_to_mat(const std::vector<std::string> &variables,
                                             const std::vector<Parameters> &param_vec)
{
  if (param_vec.size()==0)
    return arma::mat(0,0);
  
  arma::rowvec first_vector = param_vec.begin()->get_rowvec(variables);
  
  arma::mat output(param_vec.size(),first_vector.n_cols);
  output.row(0) = first_vector;
  for (size_t i=1; i<param_vec.size(); ++i)
  {
    output.row(i) = param_vec[i].get_rowvec(variables);
  }
  return output;
}

inline arma::mat vector_of_parameters_to_mat(const std::string &variable,
                                             const std::vector<Parameters> &param_vec)
{
  if (param_vec.size()==0)
    return arma::mat(0,0);
  
  arma::rowvec first_vector = param_vec.begin()->get_rowvec(variable);
  
  arma::mat output(param_vec.size(),first_vector.n_rows);
  output.row(0) = first_vector;
  for (size_t i=1; i<param_vec.size(); ++i)
  {
    output.row(i) = param_vec[i].get_rowvec(variable);
  }
  return output;
}

inline double equal_weight_gain(size_t n)
{
  return 1.0/double(n);
}

//std::vector<double> operator*(const std::vector<double> &first,
//                              const std::vector<double> &second);

//std::vector<double> operator-(const std::vector<double> &first,
//                              const std::vector<double> &second);

//std::vector<double> abs(const std::vector<double> &first);

//std::vector<double> log(const std::vector<double> &first);

//std::vector<double> operator*(const double &c,
//                              const std::vector<double> &first);

//double cess(const arma::colvec &log_weights,
//            const arma::colvec &log_incremental_weights);

// NumericMatrix get_first_element_of_list_as_numeric_matrix(const List &a_list);
// SEXP store_get_first_element_of_list_as_numeric_matrix();
//
// arma::colvec make_vector_statistic(const NumericMatrix &simulated);
// SEXP store_make_vector_statistic();

inline double inf_function(const Parameters &inputs)
{
  return arma::datum::inf;
}

inline double annealing_power(const Parameters &inputs,
                              const std::string &power_variable)
{
  return inputs[power_variable][0];
}

inline double annealing_one_minus_power(const Parameters &inputs,
                                        const std::string &power_variable)
{
  return 1.0 - inputs[power_variable][0];
}

inline arma::mat get_sigma_points(const arma::colvec &posterior_mean,
                                  const arma::mat &posterior_covariance,
                                  double w0)
{
  size_t dimension = posterior_mean.n_rows;
  arma::mat cholesky = arma::chol(posterior_covariance);
  arma::mat sigma_points = arma::mat(dimension,2*dimension+1);
  sigma_points.col(0) = posterior_mean;
  for (size_t i=1; i<dimension+1; ++i)
  {
    sigma_points.col(i) = posterior_mean + sqrt(double(dimension)/(1.0-w0))*cholesky.col(i);
  }
  for (size_t i=1; i<dimension+1; ++i)
  {
    sigma_points.col(dimension+i) = posterior_mean - sqrt(double(dimension)/(1.0-w0))*cholesky.col(i);
  }
  return(sigma_points);
}

inline arma::colvec get_unscented_weights(const arma::colvec &mean,
                                          double w0)
{
  size_t dimension = mean.n_rows;
  double weight = (1.0-w0)/double(2*dimension);
  arma::colvec unscented_weights(2*dimension+1,arma::fill::value(weight));
  unscented_weights[0] = w0;
  return(unscented_weights);
}

inline Parameters identity_transform_function(const Parameters &parameters)
{
  return parameters;
}

inline double scaled_euclidean_distance(const arma::colvec &x,
                                        const arma::colvec &y,
                                        const arma::colvec &s)
{
  return sqrt(arma::sum(arma::pow((x-y)/s,2.0)));
}

inline double scaled_Lp_distance(const arma::colvec &x,
                                 const arma::colvec &y,
                                 const arma::colvec &s,
                                 double p)
{
  return pow(arma::sum(arma::pow((x-y)/s,p)),1/p);
}

/*
inline double uniform_abc_kernel(double simulated,
                                 double observed,
                                 double tolerance)
{
  return dunif(simulated,observed-tolerance,observed+tolerance);
}

inline double uniform_abc_kernel(const arma::colvec &simulated,
                                 const arma::colvec &observed,
                                 double tolerance)
{
  return dunif(simulated,observed-tolerance,observed+tolerance);
}

inline double uniform_abc_kernel(const arma::colvec &simulated,
                                 const arma::colvec &observed,
                                 const arma::colvec &tolerances)
{
  return dunif(simulated,observed-tolerances,observed+tolerances);
}

inline double uniform_abc_kernel(const arma::mat &simulated,
                                 const arma::mat &observed,
                                 double tolerance)
{
  return dunif(simulated,observed-tolerance,observed+tolerance);
}

inline double uniform_abc_kernel(const arma::mat &simulated,
                                 const arma::mat &observed,
                                 const arma::mat &tolerances)
{
  return dunif(simulated,observed-tolerances,observed+tolerances);
}

inline double a_uniform_abc_kernel(const Parameters &simulated_data,
                                 const Parameters &observed_data,
                                 const std::string &data_name,
                                 const std::string &tolerance_name)
{
  return uniform_abc_kernel(simulated_data[data_name],observed_data[data_name],simulated_data[tolerance_name][0]);
}

inline auto make_uniform_abc_kernel(const std::string &data_name,
                                                        const std::string &tolerance_name)
{
  using namespace std::placeholders;
  return std::bind(a_uniform_abc_kernel, _1, _2, data_name, tolerance_name);
}
*/

/*
inline double evaluate_log_gaussian_random_walk_proposal(const Parameters &proposed_parameters,
 const Parameters &old_parameters,
 const Parameters &proposal_parameters)
{

}

inline Parameters simulate_gaussian_random_walk_proposal(RandomNumberGenerator &rng,
                                                         const Parameters &conditioned_on_parameters,
                                                         const Parameters &proposal_parameters)
{
  vector_parameter_iterator 
}
*/

/*
double Lp_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
                                          const double &abc_tolerance,
                                          const double &p);

double L1_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
                                          const double &abc_tolerance);
SEXP store_L1_uniform_evaluate_log_abc_kernel();

double L2_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                          const arma::colvec &observed_summary_stats,
                                          const double &abc_tolerance);
SEXP store_L2_uniform_evaluate_log_abc_kernel();

double Linf_uniform_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                            const arma::colvec &observed_summary_stats,
                                            const double &abc_tolerance);
SEXP store_Linf_uniform_evaluate_log_abc_kernel();

double gaussian_evaluate_log_abc_kernel(const arma::colvec &simulated_summary_stats,
                                        const arma::colvec &observed_summary_stats,
                                        const double &abc_tolerance);
SEXP store_gaussian_uniform_evaluate_log_abc_kernel();
 */

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

#endif
