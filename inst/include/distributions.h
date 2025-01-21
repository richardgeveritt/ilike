/**
 * @file distributions.h
 * @brief Global functions for simulation from, and evaluating the pdf/pmf of
 * distributions.
 */

#ifndef SIMULATION
#define SIMULATION

#include <RcppArmadillo.h>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <chrono>
#include <dqrng_generator.h>
#include <math.h>

// save compiler switches
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <boost/math/distributions.hpp>
// restore compiler switches
#pragma GCC diagnostic pop

#define BOOST_DISABLE_ASSERTS 1

using namespace std;

namespace ilike {

/*!
A typedef for a random number generator, using objects from dqrng.
*/
typedef dqrng::random_64bit_wrapper<dqrng::xoshiro256plus>
    RandomNumberGenerator;

/*!
A function to read the time stamp counter that works across platforms.
*/
inline size_t rdtsc() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

/*!
Generates a random number from the uniform distribution between 0 and 1.
*/
inline double runif(RandomNumberGenerator &rng) {
  boost::random::uniform_real_distribution<double> my_uniform;
  return my_uniform(rng);
}

/*!
Generates a random number from the uniform distribution between `lower` and
`upper`.
*/
inline double runif(RandomNumberGenerator &rng, double lower, double upper) {
  boost::random::uniform_real_distribution<double> my_uniform(lower, upper);
  return my_uniform(rng);
}

/*!
Generates a matrix of random numbers, where the $(i,j)$th element is drawn from
the uniform distribution between the $(i,j)$th element of `lower` and the
$(i,j)$th element of `upper`.
*/
inline arma::mat runif(RandomNumberGenerator &rng, const arma::mat &lower,
                       const arma::mat &upper) {
  arma::mat output(arma::size(lower));
  for (size_t j = 0; j < lower.n_cols; ++j) {
    for (size_t i = 0; i < lower.n_rows; ++i) {
      output.at(i, j) = runif(rng, lower.at(i, j), upper.at(i, j));
    }
  }
  return output;
}

/*!
Evaluates the log of the uniform density between 0 and 1 at `x`.
*/
inline double dunif(double x) {
  if ((x < 0.0) || (x > 1.0))
    return -arma::datum::inf;
  else {
    return 0.0;
  }
}

/*!
Evaluates the log of the uniform density between `lower` and `upper` at `x`.
*/
inline double dunif(double x, double lower, double upper) {
  if ((x < lower) || (x > upper))
    return -arma::datum::inf;
  else {
    return -log(upper - lower);
  }
}

/*!
For the $(i,j)$th element of `x`, evaluate the log of the uniform density
between the $(i,j)$th element of `lower` and the $(i,j)$th element of `upper`,
then take the sum of the result (i.e. the product of independent uniform
densities).
*/
inline double dunif(const arma::mat &x, const arma::mat &lower,
                    const arma::mat &upper) {
  double output = 0.0;
  for (size_t j = 0; j < x.n_cols; ++j) {
    for (size_t i = 0; i < x.n_rows; ++i) {
      output = output + dunif(x.at(i, j), lower.at(i, j), upper.at(i, j));
    }
  }
  return output;
}

/*!
Generates `n` random numbers from the uniform distribution between 0 and 1, and
stores them in a row vector.
*/
inline arma::rowvec multiple_runif(RandomNumberGenerator &rng, size_t n) {
  boost::random::uniform_real_distribution<double> my_uniform;
  arma::rowvec output(n);
  for (auto i = output.begin(); i != output.end(); ++i) {
    *i = my_uniform(rng);
  }
  return output;
}

/*!
Generates `n` random numbers from the uniform distribution between `lower` and
`upper`, and stores them in a row vector.
*/
inline arma::rowvec multiple_runif(RandomNumberGenerator &rng, size_t n,
                                   double lower, double upper) {
  boost::random::uniform_real_distribution<double> my_uniform(lower, upper);
  arma::rowvec output(n);
  for (auto i = output.begin(); i != output.end(); ++i) {
    *i = my_uniform(rng);
  }
  return output;
}

/*!
Generates a random number from the discrete distribution taking values [0,n),
with probabilities given by the vector `probabilities`, and `n` being taken to
be the length of this vector. The vector `probabilities` need not be normalised.
*/
inline size_t rdis(RandomNumberGenerator &rng,
                   const arma::colvec &probabilities) {
  boost::random::discrete_distribution<size_t> my_discrete(probabilities);
  return my_discrete(rng);
}

/*!
Generates `n` random numbers from the discrete distribution taking values [0,n),
with probabilities given by the vector `probabilities`, and `n` being taken to
be the length of this vector. The vector `probabilities` need not be normalised.
*/
inline arma::uvec multiple_rdis(RandomNumberGenerator &rng, size_t n,
                                const arma::colvec &probabilities) {
  boost::random::discrete_distribution<size_t> my_discrete(probabilities);
  arma::uvec output(n);
  for (auto i = output.begin(); i != output.end(); ++i) {
    *i = my_discrete(rng);
  }
  return output;
}

/*!
Generates a random number from the exponential distribution with rate `rate`.
*/
inline double rexp(RandomNumberGenerator &rng, double rate) {
  boost::random::exponential_distribution<double> my_exponential(rate);
  return my_exponential(rng);
}

/*!
Generates `n` random numbers from the exponential distribution with rate `rate`.
*/
inline arma::colvec rexp(RandomNumberGenerator &rng, size_t n, double rate) {
  boost::random::exponential_distribution<double> my_exponential(rate);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_exponential(rng);
  }
  return output;
}

/*!
Evaluates the log of the exponential density with rate `rate` at the point `x`.
*/
inline double dexp(double x, double rate) {
  if ((rate <= 0))
    return double(NAN);
  if (x < 0)
    return -arma::datum::inf;
  return log(rate) - rate * x;
}

/*!
Generates a random number from the translated (by `min`) exponential
distribution with rate `rate`.
*/
inline double rtranslatedexp(RandomNumberGenerator &rng, double rate,
                             double min) {
  boost::random::exponential_distribution<double> my_exponential(rate);
  return my_exponential(rng) + min;
}

/*!
Generates `n` random numbers from the translated (by `min`) exponential
distribution with rate `rate`.
*/
inline arma::colvec rtranslatedexp(RandomNumberGenerator &rng, size_t n,
                                   double rate, double min) {
  boost::random::exponential_distribution<double> my_exponential(rate);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_exponential(rng) + min;
  }
  return output;
}

/*!
Evaluates the log of the translated exponential (by `min`) density with rate
`rate` at the point `x`.
*/
inline double dtranslatedexp(double x, double rate, double min) {
  if ((rate <= 0))
    return double(NAN);
  if (x < min)
    return -arma::datum::inf;
  return log(rate) - rate * (x - min);
}

/*!
Generates a random number from the Poisson distribution with rate `rate`.
*/
inline int rpois(RandomNumberGenerator &rng, double rate) {
  boost::random::poisson_distribution<int> my_poisson(rate);
  return my_poisson(rng);
}

/*!
Generates `n` random numbers from the Poisson distribution with rate `rate`.
*/
inline arma::colvec rpois(RandomNumberGenerator &rng, size_t n, double rate) {
  boost::random::poisson_distribution<int> my_poisson(rate);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_poisson(rng);
  }
  return output;
}

/*!
Evaluates the log of the Poisson mass function with rate `rate` at the point
`x`.
*/
inline double dpois(double x, double rate) {
  if ((rate <= 0))
    return double(NAN);
  if (x < 0)
    return -arma::datum::inf;
  return x * log(rate) - rate - std::lgamma(x + 1);
}

/*!
Generates a random number from the standard normal distribution.
*/
inline double rnorm(RandomNumberGenerator &rng) {
  boost::random::normal_distribution<double> my_normal(0.0, 1.0);
  return my_normal(rng);
}

/*!
 Generates a matrix (of dimension given by `dimensions`) of random numbers from
 the standard normal distribution.
 */
inline arma::mat rnorm(RandomNumberGenerator &rng,
                       std::pair<size_t, size_t> dimensions) {
  boost::random::normal_distribution<double> my_normal(0.0, 1.0);
  arma::mat output(dimensions.first, dimensions.second);
  for (size_t j = 0; j < dimensions.second; ++j) {
    for (size_t i = 0; i < dimensions.first; ++i) {
      output(i, j) = my_normal(rng);
    }
  }
  return output;
}

/*!
 Generates a matrix (of dimension given by `dimensions`) of random numbers from
 the normal distribution with mean `mean` and standard deviation `sd`.
*/
inline arma::mat rnorm(RandomNumberGenerator &rng,
                       std::pair<size_t, size_t> dimensions, double mean,
                       double sd) {
  boost::random::normal_distribution<double> my_normal(mean, sd);
  arma::mat output(dimensions.first, dimensions.second);
  for (size_t j = 0; j < dimensions.second; ++j) {
    for (size_t i = 0; i < dimensions.first; ++i) {
      output(i, j) = my_normal(rng);
    }
  }
  return output;
}

/*!
Generates `n` random numbers from the standard normal distribution.
*/
inline arma::mat rnorm(RandomNumberGenerator &rng, size_t n) {
  boost::random::normal_distribution<double> my_normal(0.0, 1.0);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_normal(rng);
  }
  return output;
}

/*!
Generates `n` random numbers from the normal distribution with mean `mean` and
standard deviation `sd`.
*/
inline arma::mat rnorm(RandomNumberGenerator &rng, size_t n, double mean,
                       double sd) {
  boost::random::normal_distribution<double> my_normal(mean, sd);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_normal(rng);
  }
  return output;
}

/*!
Generates a random number from the normal distribution with mean `mean` and
standard deviation `sd`.
*/
inline double rnorm(RandomNumberGenerator &rng, double mean, double sd) {
  boost::random::normal_distribution<double> my_normal(mean, sd);
  return my_normal(rng);
}

/*!
Evaluates the log of the standard normal density at the point `x`.
*/
inline double dnorm(double x) {
  return -0.5 * log(2.0 * M_PI) - 0.5 * std::pow(x, 2.0);
}

/*!
 Evaluates the log of the normal density with mean `mean` and standard deviation
 `sd` at the point `x`.
*/
inline double dnorm(double x, double mean, double sd) {
  if (sd < 0) {
    return NAN;
  }
  if (sd == 0) {
    if (x == mean)
      return arma::datum::inf;
    else
      return -arma::datum::inf;
  }
  return -log(sd) - 0.5 * log(2.0 * M_PI) -
         0.5 * std::pow((x - mean) / sd, 2.0);
}

/*!
Evaluates the log of the normal density with mean `mean` and standard deviation
`sd` at the vector of points `x`, returning a vector. If `mean` and/or `sd` are
a `double`, the same mean/sd is used for all dimensions; if a vector the
corresponding entry in that vector is used for each dimension.
*/
inline arma::colvec dnorm(const arma::colvec &x, double mean, double sd) {
  size_t n = x.size();
  arma::colvec result(n);

  if (sd < 0) {
    for (size_t i = 0; i < n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd == 0) {
    for (size_t i = 0; i < n; ++i) {
      if (x[i] == mean) {
        result[i] = arma::datum::inf;
      } else {
        result[i] = -arma::datum::inf;
      }
    }
    return result;
  }

  for (size_t i = 0; i < n; ++i)
    result[i] = -log(sd) - 0.5 * log(2.0 * M_PI) -
                0.5 * std::pow((x[i] - mean) / sd, 2.0);
  return result;
}

/*!
Evaluates the log of the normal density with mean `mean` and standard deviation
`sd` at the vector of points `x`, returning a vector. If `mean` and/or `sd` are
a `double`, the same mean/sd is used for all dimensions; if a vector the
corresponding entry in that vector is used for each dimension.
*/
inline arma::colvec dnorm(const arma::colvec &x, const arma::colvec &mean,
                          double sd) {
  size_t n = x.size();
  arma::colvec result(n);

  if (sd < 0) {
    for (size_t i = 0; i < n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd == 0) {
    for (size_t i = 0; i < n; ++i) {
      if (x[i] == mean[i]) {
        result[i] = arma::datum::inf;
      } else {
        result[i] = -arma::datum::inf;
      }
    }
    return result;
  }

  for (size_t i = 0; i < n; ++i)
    result[i] = -log(sd) - 0.5 * log(2.0 * M_PI) -
                0.5 * std::pow((x[i] - mean[i]) / sd, 2.0);
  return result;
}

/*!
Evaluates the log of the normal density with mean `mean` and standard deviation
`sd` at the vector of points `x`, returning a vector. If `mean` and/or `sd` are
a `double`, the same mean/sd is used for all dimensions; if a vector the
corresponding entry in that vector is used for each dimension.
*/
inline arma::colvec dnorm(const arma::colvec &x, double mean,
                          const arma::colvec &sd) {
  size_t n = x.size();
  arma::colvec result(n);

  for (size_t i = 0; i < n; ++i) {
    if (sd[i] < 0) {
      result[i] = NAN;
    } else if (sd[i] == 0) {
      if (x[i] == mean) {
        result[i] = arma::datum::inf;
      } else {
        result[i] = -arma::datum::inf;
      }
    } else {
      result[i] = -log(sd[i]) - 0.5 * log(2.0 * M_PI) -
                  0.5 * std::pow((x[i] - mean) / sd[i], 2.0);
    }
  }

  return result;
}

/*!
Evaluates the log of the normal density with mean `mean` and standard deviation
`sd` at the vector of points `x`, returning a vector. If `mean` and/or `sd` are
a `double`, the same mean/sd is used for all dimensions; if a vector the
corresponding entry in that vector is used for each dimension.
*/
inline arma::colvec dnorm(const arma::colvec &x, const arma::colvec &mean,
                          const arma::colvec &sd) {
  size_t n = x.size();
  arma::colvec result(n);

  for (size_t i = 0; i < n; ++i) {
    if (sd[i] < 0) {
      result[i] = NAN;
    } else if (sd[i] == 0) {
      if (x[i] == mean[i]) {
        result[i] = arma::datum::inf;
      } else {
        result[i] = -arma::datum::inf;
      }
    } else {
      result[i] = -log(sd[i]) - 0.5 * log(2.0 * M_PI) -
                  0.5 * std::pow((x[i] - mean[i]) / sd[i], 2.0);
    }
  }

  return result;
}

/*!
Evaluates the cdf of the normal distribution with mean `mean` and standard
deviation `sd` at the point.
*/
inline double pnorm(double x, double mean, double sd) {
  return 0.5 * (1.0 + erf((x - mean) / (sd * sqrt(2.0))));
}

/*!
Generates a random number from the standard truncated normal distribution, with
bounds given by `min` and `max`.
*/
inline double rtnorm(RandomNumberGenerator &rng, double min, double max) {
  // from https://arxiv.org/abs/0907.4010

  if (min == arma::datum::inf) {
    return arma::datum::inf;
  }

  if (max == -arma::datum::inf) {
    return -arma::datum::inf;
  }

  if ((min == -arma::datum::inf) && (max == arma::datum::inf)) {
    return rnorm(rng);
  }

  if (max == arma::datum::inf) {
    if (min >= 0.0) {
      while (true) {
        double z = rnorm(rng);
        if (z >= min)
          return z;
      }
    } else {
      double alpha_star = 0.5 * (min + sqrt(min * min + 4));
      while (true) {
        double z = rtranslatedexp(rng, alpha_star, min);
        double varrho = exp(-std::pow(z - alpha_star, 2.0) / 2.0);
        if (runif(rng) <= varrho)
          return z;
      }
    }
  }

  if (min == -arma::datum::inf) {
    if (max >= 0.0) {
      while (true) {
        double z = rnorm(rng);
        if (z <= max)
          return z;
      }
    } else {
      double minus_max = -max;
      double alpha_star = 0.5 * (minus_max + sqrt(minus_max * minus_max + 4));
      while (true) {
        double z = rtranslatedexp(rng, alpha_star, minus_max);
        double varrho = exp(-std::pow(z - alpha_star, 2.0) / 2.0);
        if (runif(rng) <= varrho)
          return -z;
      }
    }
  }

  // two-sided
  while (true) {
    double z = runif(rng, min, max);
    double varrho;
    if (max < 0.0) {
      varrho = exp((max * max - z * z) / 2.0);
    } else if (0.0 < min) {
      varrho = exp((min * min - z * z) / 2.0);
    } else {
      varrho = exp(-z * z / 2.0);
    }

    if (runif(rng) <= varrho)
      return z;
  }
}

/*!
Generates a matrix (of dimension given by `dimensions`) of random numbers from
the standard truncated normal distribution, with bounds given by `min` and
`max`.
*/
inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        std::pair<size_t, size_t> dimensions, double min,
                        double max) {
  arma::mat output(dimensions.first, dimensions.second);
  for (size_t j = 0; j < dimensions.second; ++j) {
    for (size_t i = 0; i < dimensions.first; ++i) {
      output(i, j) = rtnorm(rng, min, max);
    }
  }
  return output;
}

/*!
 Generates a matrix (of dimension given by `dimensions`) of random numbers from
 the truncated normal distribution with mean `mean` and standard deviation `sd`,
 with bounds given by `min` and `max`.
*/
inline arma::mat rtnorm(RandomNumberGenerator &rng,
                        std::pair<size_t, size_t> dimensions, double mean,
                        double sd, double min, double max) {
  arma::mat output(dimensions.first, dimensions.second);
  for (size_t j = 0; j < dimensions.second; ++j) {
    for (size_t i = 0; i < dimensions.first; ++i) {
      output(i, j) = mean + rtnorm(rng, min, max) * sd;
    }
  }
  return output;
}

/*!
Generates `n` random numbers from the standard truncated normal distribution,
with bounds given by `min` and `max`.
*/
inline arma::mat rtnorm(RandomNumberGenerator &rng, size_t n, double min,
                        double max) {
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = rtnorm(rng, min, max);
  }
  return output;
}

/*!
Generates a random number from the truncated normal distribution with mean
`mean` and standard deviation `sd`, with bounds given by `min` and `max`.
*/
inline double rtnorm(RandomNumberGenerator &rng, double mean, double sd,
                     double min, double max) {
  if (sd <= 0) {
    return NAN;
  } else {
    return mean + sd * rtnorm(rng, (min - mean) / sd, (max - mean) / sd);
  }
}

/*!
Generates `n` random numbers from the truncated normal distribution with mean
`mean` and standard deviation `sd`, with bounds given by `min` and `max`.
*/
inline arma::mat rtnorm(RandomNumberGenerator &rng, size_t n, double mean,
                        double sd, double min, double max) {
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = rtnorm(rng, mean, sd, min, max);
  }
  return output;
}

/*!
Evaluates the log of the standard truncated normal density, with bounds given by
`min` and `max`, at the point `x`.
*/
inline double dtnorm(double x, double min, double max) {
  if (x < min)
    return -arma::datum::inf;
  else if (x > max)
    return -arma::datum::inf;
  else {
    return -0.5 * log(2.0 * M_PI) - 0.5 * std::pow(x, 2.0) -
           log(erf(max) - erf(min));
  }
}

/*!
Evaluates the log of the truncated normal density with mean `mean` and standard
deviation `sd`, with bounds given by `min` and `max`, at the point `x`.
*/
inline double dtnorm(double x, double mean, double sd, double min, double max) {
  if (x < min)
    return -arma::datum::inf;
  else if (x > max)
    return -arma::datum::inf;
  else {
    if (sd < 0) {
      return NAN;
    }
    if (sd == 0) {
      if (x == mean)
        return arma::datum::inf;
      else
        return -arma::datum::inf;
    }
    return -log(sd) - 0.5 * log(2.0 * M_PI) -
           0.5 * std::pow((x - mean) / sd, 2.0) -
           log(erf((max - mean) / sd) - erf((min - mean) / sd));
  }
}

/*!
Evaluates the log of the truncated normal density with mean `mean` and standard
deviation `sd`, with bounds given by `min` and `max`, at the vector of points
`x`, returning a vector.
*/
inline arma::colvec dtnorm(const arma::colvec &x, double mean, double sd,
                           double min, double max) {
  size_t n = x.size();
  arma::colvec result(n);

  if (sd < 0) {
    for (size_t i = 0; i < n; ++i)
      result[i] = NAN;
    return result;
  }
  if (sd == 0) {
    for (size_t i = 0; i < n; ++i) {
      if (x[i] < min)
        result[i] = -arma::datum::inf;
      else if (x[i] > max)
        result[i] = -arma::datum::inf;
      else {
        if (x[i] == mean) {
          result[i] = arma::datum::inf;
        } else {
          result[i] = -arma::datum::inf;
        }
      }
    }
    return result;
  }

  double erf_part = -log(erf((max - mean) / sd) - erf((min - mean) / sd));

  for (size_t i = 0; i < n; ++i) {
    if (x[i] < min)
      result[i] = -arma::datum::inf;
    else if (x[i] > max)
      result[i] = -arma::datum::inf;
    else
      result[i] = -log(sd) - 0.5 * log(2.0 * M_PI) -
                  0.5 * std::pow((x[i] - mean) / sd, 2.0) + erf_part;
  }
  return result;
}

/*!
Converts the Cholesky decomposition of a matrix to its inverse.
*/
inline arma::mat chol2inv(const arma::mat &chol) {
  arma::mat inv_chol = arma::inv(chol);
  return inv_chol * arma::trans(inv_chol);
}

/*!
Converts the Cholesky decomposition of a matrix to its log determinant.
*/
inline double chol2logdet(const arma::mat &chol) {
  return 2.0 * trace(log(chol));
}

/*!
Generates a random vector from the multivariate normal distribution with mean
`mu` and Cholesky decomposition `chol` of the covariance.
*/
inline arma::colvec rmvnorm_using_chol(RandomNumberGenerator &rng,
                                       const arma::colvec &mu,
                                       const arma::mat &chol) {
  // rmvnorm from
  // https://gallery.rcpp.org/articles/simulate-multivariate-normal/
  size_t n = 1;
  int ncols = chol.n_cols;
  arma::mat Y = rnorm(rng, std::pair<size_t, size_t>(ncols, n));
  return (mu + chol * Y);
}

/*!
Generates a random vector from the multivariate normal distribution with mean
`mu` and covariance `Sigma`.
*/
inline arma::colvec rmvnorm(RandomNumberGenerator &rng, const arma::colvec &mu,
                            const arma::mat &Sigma) {
  size_t n = 1;
  int ncols = Sigma.n_cols;
  arma::mat Y = rnorm(rng, std::pair<size_t, size_t>(ncols, n));
  return (mu + arma::chol(Sigma) * Y);
}

/*!
Generates `n` random vectors from the multivariate normal distribution with mean
`mu` and covariance `Sigma`, stored in a $d\times n$ matrix, where $d$ is the
dimension of the multivariate Gaussian (determined by the size of the inputted
mean and covariance).
*/
inline arma::mat rmvnorm(RandomNumberGenerator &rng, size_t n,
                         const arma::colvec &mu, const arma::mat &Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = rnorm(rng, std::pair<size_t, size_t>(ncols, n));
  arma::mat prod = arma::chol(Sigma) * Y;
  prod.each_col() += mu;
  return (prod);
}

/*!
Evaluates the log of the multivariate normal density with mean `mu` and
covariance `Sigma` at the vector `x`. This function does not require the
covariance as an argument: instead it takes the inverse of the covariance
`inv_Sigma` and the log of the determinant of the covariance (`log_det`). This
is designed for cases where the density needs to be evaluated at a number of
points using the same covariance; using this function avoids performing the most
expensive parts of the computation more than once.
*/
inline double dmvnorm_using_precomp(const arma::colvec &x,
                                    const arma::colvec &mu,
                                    const arma::mat &inv_Sigma,
                                    double log_det) {
  arma::colvec x_minus_mean = x - mu;
  double result =
      -((arma::size(inv_Sigma)[0] / 2.0) * log(2.0 * M_PI)) - 0.5 * log_det;
  arma::mat b = x_minus_mean.t() * inv_Sigma * x_minus_mean;
  return result - 0.5 * b(0, 0);
}

/*!
Evaluates the log of the multivariate normal density with mean `mu` and
covariance `Sigma` at the vector `x`.
*/
inline double dmvnorm(const arma::colvec &x, const arma::colvec &mu,
                      const arma::mat &Sigma) {
  double result;
  arma::colvec x_minus_mean = x - mu;
  result = -((arma::size(Sigma)[0] / 2.0) * log(2.0 * M_PI)) -
           0.5 * arma::log_det_sympd(Sigma);

  arma::mat b = x_minus_mean.t() * arma::inv_sympd(Sigma) * x_minus_mean;
  result = result - 0.5 * b(0, 0);
  return result;
}

/*!
Generates a random vector from the truncated multivariate normal distribution
with mean `mu` and Cholesky decomposition `chol` of the covariance., with bounds
given by `min` and `max`. A simple rejection sampler is used for the
implementation.
*/
inline arma::colvec rtmvnorm_using_chol(RandomNumberGenerator &rng,
                                        const arma::colvec &mu,
                                        const arma::mat &chol,
                                        const arma::colvec &min,
                                        const arma::colvec &max) {
  while (true) {
    arma::colvec z = rmvnorm_using_chol(rng, mu, chol);
    if (arma::all(z >= min) && arma::all(z <= max))
      return z;
  }
}

/*!
Generates a random vector from the truncated multivariate normal distribution
with mean `mu` and covariance `Sigma`, with bounds given by `min` and `max`. A
simple rejection sampler is used for the implementation.
*/
inline arma::colvec rtmvnorm(RandomNumberGenerator &rng, const arma::colvec &mu,
                             const arma::mat &Sigma, const arma::colvec &min,
                             const arma::colvec &max) {
  while (true) {
    arma::colvec z = rmvnorm(rng, mu, Sigma);
    if (arma::all(z >= min) && arma::all(z <= max))
      return z;
  }
}

/*!
Generates `n` random vectors from the truncated multivariate normal distribution
with mean `mu` and covariance `Sigma`, with bounds given by `min` and `max`,
stored in a $d\times n$ matrix, where $d$ is the dimension of the multivariate
Gaussian (determined by the size of the inputted mean and covariance).
*/
inline arma::mat rtmvnorm(RandomNumberGenerator &rng, size_t n,
                          const arma::colvec &mu, const arma::mat &Sigma,
                          const arma::colvec &min, const arma::colvec &max) {
  int ncols = Sigma.n_cols;
  arma::mat Y(ncols, n);
  for (size_t i = 0; i < n; ++i) {
    Y.col(i) = rtmvnorm(rng, mu, Sigma, min, max);
  }
  return Y;
}

/*
inline double dtmvnorm(const arma::colvec &x, const arma::colvec &mu,
                       const arma::mat &Sigma, const arma::colvec &min,
                       const arma::colvec &max) {
  if (arma::any(x < min) || arma::any(x > max))
    return -arma::datum::inf;
  return dmvnorm(x, mu, Sigma) - log(pmvnorm(mu, Sigma, min, max));
}
*/

/*!
Calculation of "c", as used in unbiased estimation of a Gaussian density from an
estimated mean and covariance, from Ghurye, S. G. and Olkin, I. (1969). Unbiased
Estimation of Some Multivariate Probability Densities and Related Functions. The
Annals of Mathematical Statistics 40(4), 1261–1271.
*/
inline double log_c(size_t k, size_t nu) {
  double num = -(double(k * nu) / 2.0) * log(2.0) -
               (double(k * (k - 1)) / 4.0) * log(M_PI);
  double denom = 0.0;
  for (size_t i = 0; i < k; ++i) {
    denom = denom + lgamma(0.5 * double(nu - i + 1));
  }
  return (num - denom);
}

/*!
Unbiased estimation of a Gaussian density from an estimated mean and covariance
(from `n` points), from Ghurye, S. G. and Olkin, I. (1969). Unbiased Estimation
of Some Multivariate Probability Densities and Related Functions. The Annals of
Mathematical Statistics 40(4), 1261–1271.
*/
inline double dmvnorm_estimated_params(const arma::colvec &x,
                                       const arma::colvec &estimated_mean,
                                       const arma::mat &estimated_covariance,
                                       size_t n) {
  size_t d = estimated_mean.n_rows;
  double cormac = -(double(d) / 2.0) * log(2.0 * M_PI) + log_c(d, n - 2) -
                  log_c(d, n - 1) -
                  (double(d) / 2.0) * log(1.0 - (1.0 / double(n))) -
                  (double(n - d - 2) / 2.0) *
                      arma::log_det_sympd(double(n - 1) * estimated_covariance);
  arma::colvec x_minus_mean = x - estimated_mean;
  arma::mat inner_bracket =
      double(n - 1) * estimated_covariance -
      x_minus_mean * x_minus_mean.t() / (1.0 - (1.0 / double(n)));
  double log_phi;
  if (inner_bracket.is_sympd()) {
    log_phi = arma::log_det_sympd(inner_bracket);
  } else {
    log_phi = -arma::datum::inf;
  }
  double mclaggen = (double(n - d - 3) / 2.0) * log_phi;

  return (cormac + mclaggen);
}

/*!
Generates a random number from the gamma distribution with shape `shape` and
rate `rate`.
*/
inline double rgamma(RandomNumberGenerator &rng, double shape, double rate) {
  boost::random::gamma_distribution<double> my_gamma(shape, 1.0 / rate);
  return my_gamma(rng);
}

/*!
 Generates `n` random numbers from the gamma distribution with shape `shape` and
 rate `rate`.
*/
inline arma::colvec rgamma(RandomNumberGenerator &rng, size_t n, double shape,
                           double rate) {
  boost::random::gamma_distribution<double> my_gamma(shape, 1.0 / rate);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_gamma(rng);
  }
  return output;
}

/*!
Evaluates the log of the gamma density with shape `shape` and rate `rate` at the
point `x`.
*/
inline double dgamma(double x, double shape, double rate) {
  if ((shape <= 0) || (rate <= 0))
    return double(NAN);
  if (x < 0)
    return -arma::datum::inf;
  return -boost::math::lgamma<double>(shape) + shape * log(rate) +
         (rate - 1.0) * log(x) - rate * x;
}

/*!
Generates a random number from the log-normal distribution with mean of the log
of the variable `meanlog` and standard deviation of the log of the variable
`sdlog`.
*/
inline double rlnorm(RandomNumberGenerator &rng, double meanlog, double sdlog) {
  boost::random::lognormal_distribution<double> my_gamma(meanlog, sdlog);
  return my_gamma(rng);
}

/*!
Generates `n` random numbers from the log-normal distribution with mean of the
log of the variable `meanlog` and standard deviation of the log of the variable
`sdlog`.
*/
inline arma::colvec rlnorm(RandomNumberGenerator &rng, size_t n, double meanlog,
                           double sdlog) {
  boost::random::lognormal_distribution<double> my_lognormal(meanlog, sdlog);
  arma::colvec output(n);
  for (size_t i = 0; i < n; ++i) {
    output(i) = my_lognormal(rng);
  }
  return output;
}

/*!
Evaluates the log of the log-normal density with mean of the log of the variable
`meanlog` and standard deviation of the log of the variable `sdlog` at the point
`x`.
*/
inline double dlnorm(double x, double meanlog, double sdlog) {
  if (sdlog < 0) {
    return NAN;
  }
  if (sdlog == 0) {
    if (x == meanlog)
      return arma::datum::inf;
    else
      return -arma::datum::inf;
  }
  return -log(x) - log(sdlog) - 0.5 * log(2.0 * M_PI) -
         0.5 * std::pow((log(x) - meanlog) / sdlog, 2.0);
}

/*!
Evaluates the cdf of the log-normal distribution with mean `mean` and standard
deviation `sd` at the point.
*/
inline double plnorm(double x, double meanlog, double sdlog) {
  return 0.5 * (1.0 + erf((log(x) - meanlog) / (sdlog * sqrt(2.0))));
}

/*!
Generates a random number from the multivariate log-normal distribution with
mean of the log of the variable `meanlog` and the Cholesky decomposition
`chollog` of the covariance of the log of the variable.
*/
inline arma::colvec rmvlnorm_using_chol(RandomNumberGenerator &rng,
                                        const arma::colvec &mulog,
                                        const arma::mat &chollog) {
  return exp(rmvnorm_using_chol(rng, mulog, chollog));
}

/*!
 Generates a random number from the multivariate log-normal distribution with
 mean of the log of the variable `meanlog` and covariance of the log of the
 variable `Sigmalog`.
*/
inline arma::colvec rmvlnorm(RandomNumberGenerator &rng,
                             const arma::colvec &mulog,
                             const arma::mat &Sigmalog) {
  return exp(rmvnorm(rng, mulog, Sigmalog));
}

/*!
Evaluates the log of the multivaratie log-normal density with mean of the log of
the variable `meanlog` and covariance of the log of the variable `Sigmalog` at
the point `x`. This function does not require the covariance as an argument:
instead it takes the inverse of the covariance `inv_Sigma` and the log of the
determinant of the covariance (`log_det`). This is designed for cases where the
density needs to be evaluated at a number of points using the same covariance;
using this function avoids performing the most expensive parts of the
computation more than once.
*/
inline double dmvlnorm_using_precomp(const arma::colvec &x,
                                     const arma::colvec &mulog,
                                     const arma::mat &inv_Sigmalog,
                                     double log_det) {
  arma::colvec x_minus_mean = log(x) - mulog;
  double result =
      -((arma::size(inv_Sigmalog)[0] / 2.0) * log(2.0 * M_PI)) - 0.5 * log_det;
  arma::mat b = x_minus_mean.t() * inv_Sigmalog * x_minus_mean;
  return result + sum(log(1.0 / x)) - 0.5 * b(0, 0);
}

/*!
Evaluates the log of the multivaratie log-normal density with mean of the log of
the variable `meanlog` and covariance of the log of the variable `Sigmalog` at
the point `x`.
*/
inline double dmvlnorm(const arma::colvec &x, const arma::colvec &mulog,
                       const arma::mat &Sigmalog) {
  double result;
  arma::colvec x_minus_mean = x - mulog;
  result = -((arma::size(Sigmalog)[0] / 2.0) * log(2.0 * M_PI)) -
           0.5 * arma::log_det_sympd(Sigmalog);
  arma::mat b = x_minus_mean.t() * arma::inv_sympd(Sigmalog) * x_minus_mean;
  result = result + sum(log(1.0 / x)) - 0.5 * b(0, 0);
  return result;
}

} // namespace ilike

#endif
