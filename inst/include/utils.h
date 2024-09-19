/**
 * @file utils.h
 * @brief Miscellaneous global functions.
 */

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <functional>
#include "parameters.h"
#include "distributions.h"

using boost::any_cast;

namespace ilike
{

/*!
This function implements the LogSumExp trick on the inputted vector.
*/
inline double log_sum_exp(const arma::colvec &log_weights)
{
  size_t which_max = log_weights.index_max();
  double xmax = log_weights(which_max);
  if (xmax==-arma::datum::inf)
    return -arma::datum::inf;
  
  arma::colvec log_weights_copy = log_weights;
  log_weights_copy.shed_row(which_max);
  
  double the_sum = arma::sum(arma::exp(log_weights_copy-xmax));
  return(xmax+log1p(the_sum));
}

/*!
Stratified resampling on the given vector of `log_weights`, using the given vector of uniform random numbers in `resampling_variables`.
*/
inline std::vector<size_t> stratified_resample(const arma::colvec &log_weights,
                                               const arma::colvec &resampling_variables)
{
  size_t n = log_weights.size();
  std::vector<size_t> indices(n);
  
  arma::colvec norm_log_weights = log_weights - log_sum_exp(log_weights);
  arma::colvec W = exp(norm_log_weights);
  arma::colvec cw = cumsum(W / sum(W));
  arma::colvec resampling_variables_over_n = resampling_variables/double(n);
  
  arma::colvec seq(n);
  for (size_t i=0; i<n; ++i)
  {
    seq[i] = i;
  }
  
  arma::colvec v = resampling_variables_over_n + seq/double(n);
  
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

/*!
Computes the weighted (using weight `wt`) mean of the data found in the rows of `x`. The vector `wt` needs to be normalised.
*/
inline arma::rowvec mean_wt(const arma::mat &x,
                            const arma::colvec &wt)
{
  if (wt.n_rows!=x.n_rows)
    Rcpp::stop("meanwt - matrix and weight sizes don't match.");
  arma::mat wt_mat = repmat(wt,1,x.n_cols);
  return arma::sum(wt_mat % x,0);
}

/*!
Computes the weighted (using weight `wt`) covariance of the data found in the rows of `x`. The vector `wt` needs to be normalised.
*/
inline arma::mat cov_wt(const arma::mat &x,
                        const arma::colvec &wt)
{
  arma::rowvec mean = mean_wt(x,wt);
  arma::mat x_minus_mean = x-repmat(mean,x.n_rows,1);
  arma::mat wt_mat = repmat(wt,1,x.n_cols);
  return ((wt_mat % x_minus_mean).t()*x_minus_mean)/(1 - sum(wt%wt));
}

/*!
Extracts from `param_vec` the matrices given by the inputted `variables`, converts each to a row vector, and stores them in a matrix with one per row.
*/
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

/*!
Extracts from `param_vec` the matrices given by the inputted `variable`, converts each to a row vector, and stores them in a matrix with one per row.
*/
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

/*!
Gives 1/n.
*/
inline double equal_weight_gain(size_t n)
{
  return 1.0/double(n);
}

/*!
Returns `arma::datum::inf`.
*/
inline double inf_function(const Parameters &inputs)
{
  return arma::datum::inf;
}

/*!
Extracts from `inputs` the first element of the matrix given by the index `power_variable`.
*/
inline double annealing_power(const Parameters &inputs,
                              const std::string &power_variable)
{
  return inputs[power_variable][0];
}

/*!
Extracts from `inputs` the first element of the matrix given by the index `power_variable`, then takes this away from 1.0.
*/
inline double annealing_one_minus_power(const Parameters &inputs,
                                        const std::string &power_variable)
{
  return 1.0 - inputs[power_variable][0];
}

/*!
Gives a matrix of `sigma_points` (for use in an unscented Kalman filter) given a `mean` and `covariance`.
*/
inline arma::mat get_sigma_points(const arma::colvec &mean,
                                  const arma::mat &covariance,
                                  double w0)
{
  size_t dimension = mean.n_rows;
  arma::mat cholesky = arma::chol(covariance);
  arma::mat sigma_points = arma::mat(dimension,2*dimension+1);
  sigma_points.col(0) = mean;
  for (size_t i=1; i<dimension+1; ++i)
  {
    sigma_points.col(i) = mean + sqrt(double(dimension)/(1.0-w0))*cholesky.col(i);
  }
  for (size_t i=1; i<dimension+1; ++i)
  {
    sigma_points.col(dimension+i) = mean - sqrt(double(dimension)/(1.0-w0))*cholesky.col(i);
  }
  return(sigma_points);
}

/*!
Gives a vector of weights (for use in an unscented Kalman filter) given a `mean` and the 0th weight `w_0`.
*/
inline arma::colvec get_unscented_weights(const arma::colvec &mean,
                                          double w0)
{
  size_t dimension = mean.n_rows;
  double weight = (1.0-w0)/double(2*dimension);
  arma::colvec unscented_weights(2*dimension+1,arma::fill::value(weight));
  unscented_weights[0] = w0;
  return(unscented_weights);
}

/*!
The identity function.
*/
inline Parameters identity_transform_function(const Parameters &parameters)
{
  return parameters;
}

/*!
The Euclidean distance between `x` and `y`, scaled in each dimension by `s`.
*/
inline double scaled_euclidean_distance(const arma::colvec &x,
                                        const arma::colvec &y,
                                        const arma::colvec &s)
{
  return sqrt(arma::sum(arma::pow((x-y)/s,2.0)));
}

/*!
The Lp distance between `x` and `y`, scaled in each dimension by `s`.
*/
inline double scaled_Lp_distance(const arma::colvec &x,
                                 const arma::colvec &y,
                                 const arma::colvec &s,
                                 double p)
{
  return std::pow(arma::sum(arma::pow((x-y)/s,p)),1/p);
}

/*!
Converts an `ilike::Parameters` object to an `Rcpp::List`.
*/
// [[Rcpp::export]]
inline Rcpp::List parameters_to_list(const Parameters &parameters)
{
  Rcpp::List output;
  
  for (auto it=parameters.vector_begin();it!=parameters.vector_end();++it)
  {
    output[it->first] = *it->second.first;
  }
  
  for (auto it2=parameters.any_begin();it2!=parameters.any_end();++it2)
  {
    output[it2->first] = boost::any_cast<SEXP>(*it2->second.first);
  }
  
  return output;
}

/*!
Converts an `Rcpp::List` to an `ilike::Parameters` object.
*/
// [[Rcpp::export]]
inline Parameters list_to_parameters(const Rcpp::List &list)
{
  Parameters output;
  
  if (list.length()!=0)
  {
    Rcpp::CharacterVector names = list.names();
    
    for (size_t i=0; i<names.size(); ++i)
    {
      if (Rf_isVector(list[i]))
      {
        output[Rcpp::as<std::string>(names[i])] = arma::mat(Rcpp::as<arma::vec>(Rcpp::as<NumericVector>(list[i])));
      }
      else if (Rf_isMatrix(list[i]))
      {
        output[Rcpp::as<std::string>(names[i])] = Rcpp::as<arma::mat>(list[i]);
      }
      else
      {
        output(Rcpp::as<std::string>(names[i]))= list[i];
      }
      
    }
  }
  
  return output;
}

/*!
Implements Henze-Zirkler’s test for multivariate linearity: a translation into C++ of the code found in the MVN R package.
*/
inline double hz(const arma::mat &data_in)
{
  // A C++ translation of the hz function in the MVN package in R.
  // See https://journal.r-project.org/archive/2014-2/korkmaz-goksuluk-zararsiz.pdf
  
  double tol = 1e-25;
  
  // Something I added...
  // Remove dimensions where there is no variance.
  std::vector<int> col_indices;
  for (size_t j=0; j<data_in.n_cols; ++j)
  {
    if (arma::var(data_in.col(j))>=tol)
      col_indices.push_back(j);
  }
  
  std::vector<int> row_indices;
  row_indices.reserve(data_in.n_rows);
  for (size_t i=0; i<data_in.n_rows; ++i)
  {
    row_indices.push_back(i);
  }
  
  arma::uvec col_ind = arma::conv_to<arma::uvec>::from(col_indices);
  arma::uvec row_ind = arma::conv_to<arma::uvec>::from(row_indices);
  
  arma::mat data = data_in.submat(row_ind,col_ind);
  
  // Now back to the function from MVN.
  
  size_t n = data.n_rows;
  size_t p = data.n_cols;
  
  arma::mat S = arma::cov(data,1);
  arma::rowvec column_means = arma::mean(data,0);
  arma::mat esch_row_is_column_mean(n,p);
  arma::mat dif = data;
  for (size_t j=0; j<p; ++j)
  {
    arma::vec all_col_mean_of_col_j(n);
    all_col_mean_of_col_j.fill(column_means[j]);
    dif.col(j) = data.col(j) - all_col_mean_of_col_j;
  }
  
  arma::vec Dj = arma::diagvec(dif*(arma::inv(S,arma::inv_opts::allow_approx))*dif.t()); // squared-Mahalanobis' distances
  arma::mat X;
  arma::solve(X,S,arma::eye(arma::size(S)));
  arma::mat Y = data*X*data.t();
  arma::mat Djk = - 2.0*Y.t() + arma::diagvec(Y.t())*arma::mat(1,n,arma::fill::ones) + arma::mat(n,1,arma::fill::ones)*(arma::diagvec(Y.t()).t());
  
  double b = 1.0/(sqrt(2.0)) * std::pow(((2.0*double(p) + 1.0)/4.0),(1.0/(double(p) + 4.0))) * (std::pow(double(n),(1.0/(double(p) + 4.0)))); // smoothing parameter
  
  double HZ;
  if (arma::rank(S) == p)
  {
    HZ = double(n) * (1.0/(std::pow(double(n),2.0)) * arma::accu(exp( - (std::pow(b,2.0))/2.0 * Djk)) - 2.0 *
                      ( std::pow((1.0 + (std::pow(b,2))),( - double(p)/2.0) ) ) * (1.0/double(n)) * (arma::sum(exp( - ((std::pow(b,2.0))/(2.0 * (1.0 + (std::pow(b,2.0))))) * Dj))) + ( std::pow((1.0 + (2.0 * (std::pow(b,2.0)))),( - double(p)/2.0) ) ));
  }
  else
  {
    HZ = double(n)*4.0;
  }
  
  double wb = (1.0 + std::pow(b,2.0))*(1.0 + 3.0*std::pow(b,2.0));
  
  double a = 1.0 + 2.0*std::pow(b,2.0);
  
  double mu = 1.0 - std::pow(a,(- double(p)/2.0))*(1.0 + double(p)*std::pow(b,2.0)/a + (double(p)*(double(p) + 2.0)*(std::pow(b,4.0)))/(2.0*std::pow(a,2.0))); // HZ mean
  
  double si2 = 2.0*std::pow((1.0 + 4.0*std::pow(b,2.0)),(- double(p)/2.0)) + 2.0*std::pow(a,( - double(p))) * (1.0 + (2.0*double(p)*std::pow(b,4.0))/std::pow(a,2.0) + (3.0*double(p)*(double(p) + 2.0)*std::pow(b,8.0))/(4.0*std::pow(a,4.0))) - 4.0*std::pow(wb,( - double(p)/2.0))*(1.0 + (3.0*double(p)*std::pow(b,4.0))/(2.0*wb) + (double(p)*(double(p) + 2.0)*std::pow(b,8.0))/(2.0*std::pow(wb,2.0))); //HZ variance
  
  double pmu = log(sqrt(std::pow(mu,4.0)/(si2 + std::pow(mu,2.0)))); // lognormal HZ mean
  double psi = sqrt(log((si2 + std::pow(mu,2.0))/std::pow(mu,2.0))); // lognormal HZ variance
  
  return 1.0 - plnorm(HZ,pmu,psi); // P-value associated to the HZ statistic
}

}

#endif
