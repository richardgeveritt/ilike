#ifndef SIMULATION
#define SIMULATION

//#include <random>
//#include <boost/random/binomial_distribution.hpp>
//#include <boost/random/mersenne_twister.hpp>
//#include <boost/random/normal_distribution.hpp>
#include <math.h>
//#include <stdexcept>
#include <dqrng_distribution.h>
#include <dqrng_generator.h>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/normal_distribution.hpp>

//save compiler switches
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/math/distributions.hpp>
//restore compiler switches
#pragma GCC diagnostic pop

//#include <boost/math/distributions/gamma.hpp>
//#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/gamma.hpp>

//#include <pcg_random.hpp>
//#include <xoshiro.h>

#include <RcppArmadilloForward.h>
#include <RcppCommon.h>

#define BOOST_DISABLE_ASSERTS 1

//using RandomNumberGenerator = dqrng::random_64bit_wrapper<dqrng::xoshiro256plus>;
typedef dqrng::random_64bit_wrapper<dqrng::xoshiro256plus> RandomNumberGenerator;

//typedef std::mt19937 RandomNumberGenerator;
//using Binomial = boost::random::binomial_distribution<int>;


//using Gamma = boost::math::gamma_distribution<double>;
//using RNG = dqrng::xoshiro256plus;

//inline size_t rdtsc(){
//  size_t lo,hi;
//  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
//  return ((size_t)hi << 32) | lo;
//}

#include <chrono>

inline size_t rdtsc() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

inline double runif(RandomNumberGenerator &rng)
{
  boost::random::uniform_real_distribution<double> my_uniform;
  return my_uniform(rng);
  //return dqrng::uniform01(rng);
}

inline double runif(RandomNumberGenerator &rng, double min, double max)
{
  boost::random::uniform_real_distribution<double> my_uniform(min, max);
  return my_uniform(rng);
}

inline arma::mat runif(RandomNumberGenerator &rng,
                       const arma::mat &lower,
                       const arma::mat &upper)
{
  arma::mat output(arma::size(lower));
  for (size_t j=0; j<lower.n_cols; ++j)
  {
    for (size_t i=0; i<lower.n_rows; ++i)
    {
      output.at(i,j) = runif(rng,
                             lower.at(i,j),
                             upper.at(i,j));
    }
  }
  return output;
}

inline double dunif(double point,
                    double lower_bound,
                    double upper_bound)
{
  if ( (point<lower_bound) || (point>upper_bound) )
    return -arma::datum::inf;
  else
  {
    return -log(upper_bound-lower_bound);
  }
}

inline double dunif(const arma::mat &point,
                    const arma::mat &lower_bounds,
                    const arma::mat &upper_bounds)
{
  double output = 0.0;
  for (size_t j=0; j<point.n_cols; ++j)
  {
    for (size_t i=0; i<point.n_rows; ++i)
    {
      output = output + dunif(point.at(i,j),
                              lower_bounds.at(i,j),
                              upper_bounds.at(i,j));
    }
  }
  return output;
}

inline arma::rowvec multiple_runif(RandomNumberGenerator &rng,
                                   size_t n)
{
  boost::random::uniform_real_distribution<double> my_uniform;
  arma::rowvec output(n);
  for (auto i=output.begin();
       i!=output.end();
       ++i)
  {
    *i = my_uniform(rng);
  }
  return output;
}

inline arma::rowvec multiple_runif(RandomNumberGenerator &rng,
                                   size_t n,
                                   double min,
                                   double max)
{
  boost::random::uniform_real_distribution<double> my_uniform(min, max);
  arma::rowvec output(n);
  for (auto i=output.begin();
       i!=output.end();
       ++i)
  {
    *i = my_uniform(rng);
  }
  return output;
}

inline size_t rdis(RandomNumberGenerator &rng,
                   const arma::colvec &probabilities)
{
  boost::random::discrete_distribution<size_t> my_discrete(probabilities);
  return my_discrete(rng);
}

inline arma::uvec multiple_rdis(RandomNumberGenerator &rng,
                                size_t n,
                                const arma::colvec &probabilities)
{
  boost::random::discrete_distribution<size_t> my_discrete(probabilities);
  arma::uvec output(n);
  for (auto i=output.begin();
       i!=output.end();
       ++i)
  {
    *i = my_discrete(rng);
  }
  return output;
}

inline double rexp(RandomNumberGenerator &rng, double rate)
{
  boost::random::exponential_distribution<double> my_exponential(rate);
  return my_exponential(rng);
}

inline arma::colvec rexp(RandomNumberGenerator &rng, size_t n, double rate)
{
  boost::random::exponential_distribution<double> my_exponential(rate);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_exponential(rng);
  }
  return output;
}

inline double dexp(double x, double rate)
{
  if ( (rate<=0) )
    return double(NAN);
  if (x<0)
    return -arma::datum::inf;
  return log(rate) - rate*x;
}

inline double rtranslatedexp(RandomNumberGenerator &rng, double rate, double min)
{
  boost::random::exponential_distribution<double> my_exponential(rate);
  return my_exponential(rng)+min;
}

inline arma::colvec rtranslatedexp(RandomNumberGenerator &rng, size_t n, double rate, double min)
{
  boost::random::exponential_distribution<double> my_exponential(rate);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_exponential(rng)+min;
  }
  return output;
}

inline double dtranslatedexp(double x, double rate, double min)
{
  if ( (rate<=0) )
    return double(NAN);
  if (x<min)
    return -arma::datum::inf;
  return log(rate) - rate*(x-min);
}

inline double rnorm(RandomNumberGenerator &rng)
{
  boost::random::normal_distribution<double> my_normal(0.0,1.0);
  return my_normal(rng);
}

inline arma::mat rnorm(RandomNumberGenerator &rng,
                       std::pair<size_t,size_t> dimensions)
{
  boost::random::normal_distribution<double> my_normal(0.0,1.0);
  arma::mat output(dimensions.first,dimensions.second);
  for (size_t j=0; j<dimensions.second; ++j)
  {
    for (size_t i=0; i<dimensions.first; ++i)
    {
      output(i,j) = my_normal(rng);
    }
  }
  return output;
}

inline arma::mat rnorm(RandomNumberGenerator &rng,
                       std::pair<size_t,size_t> dimensions,
                       double mean,
                       double sd)
{
  boost::random::normal_distribution<double> my_normal(mean,sd);
  arma::mat output(dimensions.first,dimensions.second);
  for (size_t j=0; j<dimensions.second; ++j)
  {
    for (size_t i=0; i<dimensions.first; ++i)
    {
      output(i,j) = my_normal(rng);
    }
  }
  return output;
}

inline arma::mat rnorm(RandomNumberGenerator &rng,
                       size_t n)
{
  boost::random::normal_distribution<double> my_normal(0.0,1.0);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_normal(rng);
  }
  return output;
}

inline arma::mat rnorm(RandomNumberGenerator &rng,
                       size_t n,
                       double mean,
                       double sd)
{
  boost::random::normal_distribution<double> my_normal(mean,sd);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_normal(rng);
  }
  return output;
}

inline double rnorm(RandomNumberGenerator &rng, double mean, double sd)
{
  boost::random::normal_distribution<double> my_normal(mean, sd);
  return my_normal(rng);
}

inline double dnorm(double x)
{
  return - 0.5*log(2.0*M_PI) - 0.5*pow(x,2.0);
}

inline double dnorm(double x, double mean, double sd)
{
  if (sd<0)
  {
    return NAN;
  }
  if (sd==0)
  {
    if (x==mean)
      return arma::datum::inf;
    else
      return -arma::datum::inf;
  }
  return -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x-mean)/sd,2.0);
}

inline arma::colvec dnorm(const arma::colvec &x, double mean, double sd)
{
  size_t n = x.size();
  arma::colvec result(n);

  if (sd<0)
  {
    for (size_t i = 0; i<n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd==0)
  {
    for (size_t i = 0; i<n; ++i)
    {
      if (x[i]==mean)
      {
        result[i] = arma::datum::inf;
      }
      else
      {
        result[i] = -arma::datum::inf;
      }
    }
    return result;
  }

  for (size_t i = 0; i<n; ++i)
    result[i] = -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x[i]-mean)/sd,2.0);
  return result;

}

inline double rtnorm(RandomNumberGenerator &rng,
                     double min,
                     double max)
{
  // from https://arxiv.org/abs/0907.4010
  
  if (min==arma::datum::inf)
  {
    return arma::datum::inf;
  }
  
  if (max==-arma::datum::inf)
  {
    return -arma::datum::inf;
  }
    
  if ( (min==-arma::datum::inf) && (max==arma::datum::inf) )
  {
    return rnorm(rng);
  }
  
  if (max==arma::datum::inf)
  {
    if (min>=0.0)
    {
      while (true)
      {
        double z = rnorm(rng);
        if (z >= min)
          return z;
      }
    }
    else
    {
      double alpha_star = 0.5*(min+sqrt(min*min+4));
      while (true)
      {
        double z = rtranslatedexp(rng,alpha_star,min);
        double varrho = exp(-pow(z-alpha_star,2.0)/2.0);
        if (runif(rng)<=varrho)
          return z;
      }
    }
  }
  
  if (min==-arma::datum::inf)
  {
    if (max>=0.0)
    {
      while (true)
      {
        double z = rnorm(rng);
        if (z <= max)
          return z;
      }
    }
    else
    {
      double minus_max = -max;
      double alpha_star = 0.5*(minus_max+sqrt(minus_max*minus_max+4));
      while (true)
      {
        double z = rtranslatedexp(rng,alpha_star,minus_max);
        double varrho = exp(-pow(z-alpha_star,2.0)/2.0);
        if (runif(rng)<=varrho)
          return -z;
      }
    }
  }
  
  // two-sided
  while (true)
  {
    double z = runif(rng,min,max);
    double varrho;
    if (max<0.0)
    {
      varrho = exp((max*max-z*z)/2.0);
    }
    else if (0.0<min)
    {
      varrho = exp((min*min-z*z)/2.0);
    }
    else
    {
      varrho = exp(-z*z/2.0);
    }
    
    if (runif(rng)<=varrho)
      return z;
  }
}

inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        std::pair<size_t,size_t> dimensions,
                        double min,
                        double max)
{
  arma::mat output(dimensions.first,dimensions.second);
  for (size_t j=0; j<dimensions.second; ++j)
  {
    for (size_t i=0; i<dimensions.first; ++i)
    {
      output(i,j) = rtnorm(rng,min,max);
    }
  }
  return output;
}

inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        std::pair<size_t,size_t> dimensions,
                        double mean,
                        double sd,
                        double min,
                        double max)
{
  arma::mat output(dimensions.first,dimensions.second);
  for (size_t j=0; j<dimensions.second; ++j)
  {
    for (size_t i=0; i<dimensions.first; ++i)
    {
      output(i,j) = mean + rtnorm(rng,min,max)*sd;
    }
  }
  return output;
}

inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        size_t n,
                        double min,
                        double max)
{
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = rtnorm(rng,min,max);
  }
  return output;
}

inline double rtnorm(RandomNumberGenerator &rng,
                     double mean,
                     double sd,
                     double min,
                     double max)
{
  if (sd<=0)
  {
    return NAN;
  }
  else
  {
    return mean + sd*rtnorm(rng,(min-mean)/sd,(max-mean)/sd);
  }
}

inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        size_t n,
                        double mean,
                        double sd,
                        double min,
                        double max)
{
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = rtnorm(rng,mean,sd,min,max);
  }
  return output;
}

inline double dtnorm(double x,
                     double min,
                     double max)
{
  if (x<min)
    return -arma::datum::inf;
  else if (x>max)
    return -arma::datum::inf;
  else
  {
    return - 0.5*log(2.0*M_PI) - 0.5*pow(x,2.0) - log(erf(max) - erf(min));
  }
}

inline double dtnorm(double x,
                     double mean,
                     double sd,
                     double min,
                     double max)
{
  if (x<min)
    return -arma::datum::inf;
  else if (x>max)
    return -arma::datum::inf;
  else
  {
    if (sd<0)
    {
      return NAN;
    }
    if (sd==0)
    {
      if (x==mean)
        return arma::datum::inf;
      else
        return -arma::datum::inf;
    }
    return -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x-mean)/sd,2.0) - log(erf((max-mean)/sd) - erf((min-mean)/sd));
  }
}

inline arma::colvec dtnorm(const arma::colvec &x,
                           double mean,
                           double sd,
                           double min,
                           double max)
{
  size_t n = x.size();
  arma::colvec result(n);
  
  if (sd<0)
  {
    for (size_t i = 0; i<n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd==0)
  {
    for (size_t i = 0; i<n; ++i)
    {
      if (x[i]<min)
        result[i] = -arma::datum::inf;
      else if (x[i]>max)
        result[i] = -arma::datum::inf;
      else
      {
        if (x[i]==mean)
        {
          result[i] = arma::datum::inf;
        }
        else
        {
          result[i] = -arma::datum::inf;
        }
      }
    }
    return result;
  }
  
  double erf_part =  - log(erf((max-mean)/sd) - erf((min-mean)/sd));
  
  for (size_t i = 0; i<n; ++i)
  {
    if (x[i]<min)
      result[i] = -arma::datum::inf;
    else if (x[i]>max)
      result[i] = -arma::datum::inf;
    else
      result[i] = -log(sd) - 0.5*log(2.0*M_PI) - 0.5*pow((x[i]-mean)/sd,2.0) + erf_part;
  }
  return result;
  
}

inline arma::mat chol2inv(const arma::mat &chol)
{
  arma::mat inv_chol = arma::inv(chol);
  return inv_chol*arma::trans(inv_chol);
}

inline double chol2logdet(const arma::mat &chol)
{
  return 2.0*trace(log(chol));
}

// rmvnorm from https://gallery.rcpp.org/articles/simulate-multivariate-normal/

inline arma::colvec rmvnorm_using_chol(RandomNumberGenerator &rng,
                                       const arma::colvec &mu,
                                       const arma::mat &chol)
{
  //size_t n = 1;
  //int ncols = sigma.n_cols;
  //arma::mat Y = rnorm(rng,
  //                    std::pair<size_t,size_t>(ncols,n));
  //if (sigma_is_chol)
  //  return (arma::repmat(mu, n, 1) + sigma*Y).as_col();
  //else
  //  return (arma::repmat(mu, n, 1) + arma::chol(sigma)*Y).as_col();
  
  size_t n = 1;
  int ncols = chol.n_cols;
  arma::mat Y = rnorm(rng,
                      std::pair<size_t,size_t>(ncols,n));
  return (mu + chol*Y);
}

inline arma::colvec rmvnorm(RandomNumberGenerator &rng,
                            const arma::colvec &mu,
                            const arma::mat &Sigma)
{
  //size_t n = 1;
  //int ncols = sigma.n_cols;
  //arma::mat Y = rnorm(rng,
  //                    std::pair<size_t,size_t>(ncols,n));
  //if (sigma_is_chol)
  //  return (arma::repmat(mu, n, 1) + sigma*Y).as_col();
  //else
  //  return (arma::repmat(mu, n, 1) + arma::chol(sigma)*Y).as_col();
  
  size_t n = 1;
  int ncols = Sigma.n_cols;
  arma::mat Y = rnorm(rng,
                      std::pair<size_t,size_t>(ncols,n));
  return (mu + arma::chol(Sigma)*Y);
}

inline arma::mat rmvnorm(RandomNumberGenerator &rng,
                         size_t n,
                         const arma::colvec &mu,
                         const arma::mat &Sigma)
{
  //size_t n = 1;
  //int ncols = sigma.n_cols;
  //arma::mat Y = rnorm(rng,
  //                    std::pair<size_t,size_t>(ncols,n));
  //if (sigma_is_chol)
  //  return (arma::repmat(mu, n, 1) + sigma*Y).as_col();
  //else
  //  return (arma::repmat(mu, n, 1) + arma::chol(sigma)*Y).as_col();
  
  int ncols = Sigma.n_cols;
  arma::mat Y = rnorm(rng,
                      std::pair<size_t,size_t>(ncols,n));
  arma::mat prod = arma::chol(Sigma)*Y;
  prod.each_col() += mu;
  return (prod);
}

inline double dmvnorm_using_precomp(const arma::colvec &x,
                                    const arma::colvec &mu,
                                    const arma::mat &inv_Sigma,
                                    double log_det)
{
  arma::colvec x_minus_mean = x-mu;
  double result = -((arma::size(inv_Sigma)[0]/2.0) * log(2.0*M_PI)) - 0.5*log_det;
  //double thing = double(x_minus_mean.t()*arma::inv_sympd(c)*x_minus_mean(0,0));
  arma::mat b = x_minus_mean.t()*inv_Sigma*x_minus_mean;
  return result - 0.5*b(0,0);
}

inline double dmvnorm(const arma::colvec &x,
                      const arma::colvec &mu,
                      const arma::mat &Sigma)
{
  double result;
  arma::colvec x_minus_mean = x-mu;
  //try
  //{
  result = -((arma::size(Sigma)[0]/2.0) * log(2.0*M_PI)) - 0.5*arma::log_det_sympd(Sigma);
  //double thing = double(x_minus_mean.t()*arma::inv_sympd(c)*x_minus_mean(0,0));
  arma::mat b = x_minus_mean.t()*arma::inv_sympd(Sigma)*x_minus_mean;
  result = result - 0.5*b(0,0);
  //}
  //catch (std::exception)
  //{
  //  Rcpp::stop("mvnormal_logpdf - covariance is not positive definite.");
  //  //Rcpp::stop("mvnormal_logpdf - covariance is not positive definite.");
  //}
  return result;
}

// rmvnorm from https://gallery.rcpp.org/articles/simulate-multivariate-normal/

inline arma::colvec rtmvnorm_using_chol(RandomNumberGenerator &rng,
                                        const arma::colvec &mu,
                                        const arma::mat &chol,
                                        const arma::colvec &min,
                                        const arma::colvec &max)
{
  while (true)
  {
    arma::colvec z = rmvnorm_using_chol(rng,mu,chol);
    if ( arma::all(z >= min) && arma::all(z <= max) )
      return z;
  }
}

inline arma::colvec rtmvnorm(RandomNumberGenerator &rng,
                             const arma::colvec &mu,
                             const arma::mat &Sigma,
                             const arma::colvec &min,
                             const arma::colvec &max)
{
  while (true)
  {
    arma::colvec z = rmvnorm(rng,mu,Sigma);
    if ( arma::all(z >= min) && arma::all(z <= max) )
      return z;
  }
}

inline arma::mat rtmvnorm(RandomNumberGenerator &rng,
                          size_t n,
                          const arma::colvec &mu,
                          const arma::mat &Sigma,
                          const arma::colvec &min,
                          const arma::colvec &max)
{
  int ncols = Sigma.n_cols;
  arma::mat Y(ncols,n);
  for (size_t i=0; i<n; ++i)
  {
    Y.col(i) = rtmvnorm(rng,mu,Sigma,min,max);
  }
  return Y;
}

/*
// dmvnorm from https://gallery.rcpp.org/articles/dmvnorm_arma/
static double const log2pi = std::log(2.0 * M_PI);

inline void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

inline arma::vec dmvnorm(const arma::rowvec &x,
                  const arma::rowvec &mean,
                  const arma::mat &sigma)
{
  using arma::uword;
  uword const n = x.n_rows,
  xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
  constants = -(double)xdim/2.0 * log2pi,
  other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }
  
  return out;
}

inline arma::vec dmvnorm(const arma::mat &x,
                  const arma::rowvec &mean,
                  const arma::mat &sigma) {
  using arma::uword;
  uword const n = x.n_rows,
  xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
  constants = -(double)xdim/2.0 * log2pi,
  other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }
  
  return out;
}
*/


inline double log_c(size_t k, size_t nu)
{
  double num = -(double(k*nu)/2.0)*log(2.0) - (double(k*(k-1))/4.0)*log(M_PI);
  double denom = 0.0;
  for (size_t i=0; i<k; ++i)
  {
    denom = denom + lgamma(0.5*double(nu-i+1));
  }
  return(num - denom);
}

inline double dmvnorm_estimated_params(const arma::colvec &x,
                                       const arma::colvec &estimated_mean,
                                       const arma::mat &estimated_covariance,
                                       size_t n)
{
  size_t d = estimated_mean.n_rows;
  //try
  //{
  double cormac = -(double(d)/2.0)*log(2.0*M_PI) + log_c(d,n-2) - log_c(d,n-1) - (double(d)/2.0)*log(1.0-(1.0/double(n))) - (double(n-d-2)/2.0)*arma::log_det_sympd(double(n-1)*estimated_covariance);
  arma::colvec x_minus_mean = x-estimated_mean;
  arma::mat inner_bracket = double(n-1)*estimated_covariance - x_minus_mean*x_minus_mean.t()/(1.0-(1.0/double(n)));
  double log_phi;
  if (inner_bracket.is_sympd())
  {
    log_phi = arma::log_det_sympd(inner_bracket);
  }
  else
  {
    log_phi = -arma::datum::inf;
  }
  double mclaggen = (double(n-d-3)/2.0)*log_phi;
  
  return(cormac + mclaggen);
  //}
  //catch(std::exception)
  //{
  //  Rcpp::stop("mvnormal_logpdf_unbiased_with_estimated_params - it might be that n<=d-3.");
  //}

}

inline double rgamma(RandomNumberGenerator &rng, double shape, double rate)
{
  boost::random::gamma_distribution<double> my_gamma(shape, 1.0/rate);
  return my_gamma(rng);
}

inline arma::colvec rgamma(RandomNumberGenerator &rng, size_t n, double shape, double rate)
{
  boost::random::gamma_distribution<double> my_gamma(shape, 1.0/rate);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_gamma(rng);
  }
  return output;
}

inline double dgamma(double x, double shape, double rate)
{
  if ( (shape<=0) || (rate<=0) )
    return double(NAN);
  if (x<0)
    return -arma::datum::inf;
  return -boost::math::lgamma<double>(shape) + shape*log(rate) + (rate-1.0)*log(x) - rate*x;
}

inline double rlnorm(RandomNumberGenerator &rng, double meanlog, double sdlog)
{
  boost::random::lognormal_distribution<double> my_gamma(meanlog, sdlog);
  return my_gamma(rng);
}

inline arma::colvec rlnorm(RandomNumberGenerator &rng, size_t n, double meanlog, double sdlog)
{
  boost::random::lognormal_distribution<double> my_lognormal(meanlog, sdlog);
  arma::colvec output(n);
  for (size_t i=0; i<n; ++i)
  {
    output(i) = my_lognormal(rng);
  }
  return output;
}

inline double dlnorm(double x, double meanlog, double sdlog)
{
  if (sdlog<0)
  {
    return NAN;
  }
  if (sdlog==0)
  {
    if (x==meanlog)
      return arma::datum::inf;
    else
      return -arma::datum::inf;
  }
  return - log(x) - log(sdlog) - 0.5*log(2.0*M_PI) - 0.5*pow((log(x)-meanlog)/sdlog,2.0);
}

inline arma::colvec rmvlnorm_using_chol(RandomNumberGenerator &rng,
                             const arma::colvec &mulog,
                             const arma::mat &chollog)
{
  //size_t n = 1;
  //int ncols = sigma.n_cols;
  //arma::mat Y = rnorm(rng,
  //                    std::pair<size_t,size_t>(ncols,n));
  //if (sigma_is_chol)
  //  return (arma::repmat(mu, n, 1) + sigma*Y).as_col();
  //else
  //  return (arma::repmat(mu, n, 1) + arma::chol(sigma)*Y).as_col();
  
  return exp(rmvnorm_using_chol(rng,mulog,chollog));
}

inline arma::colvec rmvlnorm(RandomNumberGenerator &rng,
                             const arma::colvec &mulog,
                             const arma::mat &Sigmalog)
{
  //size_t n = 1;
  //int ncols = sigma.n_cols;
  //arma::mat Y = rnorm(rng,
  //                    std::pair<size_t,size_t>(ncols,n));
  //if (sigma_is_chol)
  //  return (arma::repmat(mu, n, 1) + sigma*Y).as_col();
  //else
  //  return (arma::repmat(mu, n, 1) + arma::chol(sigma)*Y).as_col();
  
  return exp(rmvnorm(rng,mulog,Sigmalog));
}

inline double dmvlnorm_using_precomp(const arma::colvec &x,
                                     const arma::colvec &mulog,
                                     const arma::mat &inv_Sigmalog,
                                     double log_det)
{
  arma::colvec x_minus_mean = log(x)-mulog;
  double result = -((arma::size(inv_Sigmalog)[0]/2.0) * log(2.0*M_PI)) - 0.5*log_det;
  //double thing = double(x_minus_mean.t()*arma::inv_sympd(c)*x_minus_mean(0,0));
  arma::mat b = x_minus_mean.t()*inv_Sigmalog*x_minus_mean;
  //arma::colvec ones;
  //ones.ones();
  return result + sum(log(1.0/x)) - 0.5*b(0,0);
}

inline double dmvlnorm(const arma::colvec &x,
                       const arma::colvec &mulog,
                       const arma::mat &Sigmalog)
{
  double result;
  arma::colvec x_minus_mean = x-mulog;
  //try
  //{
  result = -((arma::size(Sigmalog)[0]/2.0) * log(2.0*M_PI)) - 0.5*arma::log_det_sympd(Sigmalog);
  //double thing = double(x_minus_mean.t()*arma::inv_sympd(c)*x_minus_mean(0,0));
  arma::mat b = x_minus_mean.t()*arma::inv_sympd(Sigmalog)*x_minus_mean;
  result = result + sum(log(1.0/x)) - 0.5*b(0,0);
  //}
  //catch (std::exception)
  //{
  //  Rcpp::stop("mvnormal_logpdf - covariance is not positive definite.");
  //  //Rcpp::stop("mvnormal_logpdf - covariance is not positive definite.");
  //}
  return result;
}

// namespace Rcpp {
//   SEXP wrap(const RNG& rng);
//   RNG& as(SEXP ptr_RNG);
// }

//RCPP_EXPOSED_CLASS(pcg64);

//#include <Rcpp.h>

// namespace Rcpp {
//
//   inline SEXP wrap(RNG& rng) {
//     return Rcpp::XPtr<RNG>(new RNG(&rng));
//   }
//
//   inline RNG& as(SEXP ptr_RNG) {
//     Rcpp::XPtr<RNG> ptr(ptr_RNG);
//     RNG& rng = *ptr;
//     return rng;
//   }
// }

// #include <array>
// #include <mystdint.h>
// #include <functional>
// #include <algorithm>
// #include <type_traits>
//
// //template<size_t N, int_fast8_t A, int_fast8_t B, int_fast8_t C>
// class MyRNG {
//
// public:
//
//   //using result_type = uint64_t;
//
// private:
//
//   // int_fast8_t A;
//   // int_fast8_t B;
//   // int_fast8_t C;
//   //
//   // std::array<uint64_t, 4> state;
//   //
//   // struct SplitMix {
//   //   SplitMix(const uint64_t& k) : state(k) {}
//   //
//   //   uint64_t operator() () {
//   //     uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
//   //     z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
//   //     z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
//   //     return z ^ (z >> 31);
//   //   }
//   //
//   // private:
//   //   uint64_t state;
//   // };
//   //
//   // uint64_t rotl(const uint64_t x, int k);
//   //
//   // uint64_t next();
//
// public:
//   // inline static constexpr uint64_t min() {return 0.0;};
//   // inline static constexpr uint64_t max() {return UINT64_MAX;};
//
//   MyRNG();
//   virtual ~MyRNG();
//   //xoshiro256plus(uint64_t _seed = 0x85c6ea9eb065ebeeULL);
//
//   // void seed(std::function<uint64_t(void)> rng);
//   //
//   // void seed(uint64_t _seed);
//   //
//   // uint64_t operator() ();
//   //
//   // void jump();
//   // void jump(uint64_t n);
//   // void long_jump();
//   // void long_jump(uint64_t n);
// };
//
// RCPP_EXPOSED_WRAP(MyRNG);
//
// #include <Rcpp.h>
//
//
// // inline uint64_t xoshiro256plus::rotl(const uint64_t x, int k) {
// //   return (x << k) | (x >> (64 - k));
// // }
// //
// // inline uint64_t xoshiro256plus::next() {
// //   size_t N_ = 4;
// //   const uint64_t result = state[0] + state[N_ - 1];
// //
// //   const uint64_t t = state[1] << A;
// //
// //   state[2] ^= state[0];
// //   state[3] ^= state[1];
// //   state[1] ^= state[2];
// //   state[0] ^= state[3];
// //
// //   state[2] ^= t;
// //
// //   state[3] = rotl(state[3], B);
// //
// //   return result;
// // }
//
// inline MyRNG::MyRNG()
// {
//   //seed(_seed);
//   //this->A = 17;
//   //this->B = 45;
//   //this->C = 0;
// }
//
// inline MyRNG::~MyRNG()
// {
//   //seed(_seed);
//   //this->A = 17;
//   //this->B = 45;
//   //this->C = 0;
// }
//
// // inline xoshiro256plus::xoshiro256plus(uint64_t _seed) {
// //   //seed(_seed);
// //   //this->A = 17;
// //   //this->B = 45;
// //   //this->C = 0;
// // }
//
// // inline void xoshiro256plus::seed(std::function<uint64_t(void)> rng) {
// //   std::generate(state.begin(), state.end(), rng);
// // }
// //
// // inline void xoshiro256plus::seed(uint64_t _seed) {
// //   seed(SplitMix(_seed));
// // }
// //
// // inline uint64_t xoshiro256plus::operator() () {
// //   return next();
// // }
// //
// // inline void xoshiro256plus::jump(uint64_t n) {
// //   for( ; n > 0; --n) jump();
// // }
// //
// // inline void xoshiro256plus::long_jump(uint64_t n) {
// //   for( ; n > 0; --n) long_jump();
// // }
// //
// // /* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
// //  numbers. We suggest to use its upper bits for floating-point
// //  generation, as it is slightly faster than xoshiro256**. It passes all
// //  tests we are aware of except for the lowest three bits, which might
// //  fail linearity tests (and just those), so if low linear complexity is
// //  not considered an issue (as it is usually the case) it can be used to
// //  generate 64-bit outputs, too.
// //  We suggest to use a sign test to extract a random Boolean value, and
// //  right shifts to extract subsets of bits.
// //  The state must be seeded so that it is not everywhere zero. If you have
// //  a 64-bit seed, we suggest to seed a splitmix64 generator and use its
// //  output to fill s. */
// //
// // /* This is the jump function for the generator. It is equivalent
// //  to 2^128 calls to next(); it can be used to generate 2^128
// //  non-overlapping subsequences for parallel computations. */
// // inline void xoshiro256plus::jump() {
// //   static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
// //
// //   uint64_t s0 = 0;
// //   uint64_t s1 = 0;
// //   uint64_t s2 = 0;
// //   uint64_t s3 = 0;
// //   for(unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
// //     for(unsigned int b = 0; b < 64; b++) {
// //       if (JUMP[i] & uint64_t(1) << b) {
// //         s0 ^= state[0];
// //         s1 ^= state[1];
// //         s2 ^= state[2];
// //         s3 ^= state[3];
// //       }
// //       operator()();
// //     }
// //
// //     state[0] = s0;
// //   state[1] = s1;
// //   state[2] = s2;
// //   state[3] = s3;
// // }
// //
// // /* This is the long-jump function for the generator. It is equivalent to
// //  2^192 calls to next(); it can be used to generate 2^64 starting points,
// //  from each of which jump() will generate 2^64 non-overlapping
// //  subsequences for parallel distributed computations. */
// // inline void xoshiro256plus::long_jump(void) {
// //   static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };
// //
// //   uint64_t s0 = 0;
// //   uint64_t s1 = 0;
// //   uint64_t s2 = 0;
// //   uint64_t s3 = 0;
// //   for(unsigned int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
// //     for(unsigned int b = 0; b < 64; b++) {
// //       if (LONG_JUMP[i] & uint64_t(1) << b) {
// //         s0 ^= state[0];
// //         s1 ^= state[1];
// //         s2 ^= state[2];
// //         s3 ^= state[3];
// //       }
// //       operator()();
// //     }
// //
// //     state[0] = s0;
// //   state[1] = s1;
// //   state[2] = s2;
// //   state[3] = s3;
// // }

#endif
