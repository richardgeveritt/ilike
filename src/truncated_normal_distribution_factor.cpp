#include "truncated_normal_distribution_factor.h"

namespace ilike
{

TruncatedNormalDistributionFactor::TruncatedNormalDistributionFactor()
:DistributionFactor(), mean(0.0), sd(1.0), lower(-arma::datum::inf), upper(arma::datum::inf)
{
}

TruncatedNormalDistributionFactor::TruncatedNormalDistributionFactor(const std::string &variable_in, double mean_in, double sd_in,
                                                                     double lower_in, double upper_in)
:DistributionFactor(), variable(variable_in), mean(mean_in), sd(sd_in), lower(lower_in), upper(upper_in)
{
  if (!(sd_in > 0.0) || !(lower_in < upper_in))
    Rcpp::stop("ilike::tnorm: the standard deviation must be positive and lower < upper.");
}

TruncatedNormalDistributionFactor::~TruncatedNormalDistributionFactor()
{
}

TruncatedNormalDistributionFactor::TruncatedNormalDistributionFactor(const TruncatedNormalDistributionFactor &another)
:DistributionFactor(another)
{
  this->make_copy(another);
}

void TruncatedNormalDistributionFactor::operator=(const TruncatedNormalDistributionFactor &another)
{
  if (this == &another)
    return;
  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor* TruncatedNormalDistributionFactor::duplicate() const
{
  return new TruncatedNormalDistributionFactor(*this);
}

DistributionFactor* TruncatedNormalDistributionFactor::distribution_factor_duplicate() const
{
  return new TruncatedNormalDistributionFactor(*this);
}

void TruncatedNormalDistributionFactor::make_copy(const TruncatedNormalDistributionFactor &another)
{
  this->variable = another.variable;
  this->mean = another.mean;
  this->sd = another.sd;
  this->lower = another.lower;
  this->upper = another.upper;
}

double TruncatedNormalDistributionFactor::distribution_evaluate(const Parameters &input) const
{
  arma::colvec x = input.get_colvec(this->variable);
  double output = 0.0;
  for (arma::uword k = 0; k < x.n_elem; ++k)
    output = output + log_truncated_normal(x[k], this->mean, this->sd, this->lower, this->upper);
  return output;
}

arma::mat TruncatedNormalDistributionFactor::distribution_evaluate_gradient(const std::string &variable_in,
                                                                            const Parameters &input) const
{
  // with respect to x (the arguments are numbers): -(x - mean)/sd^2 inside
  // the limits (the density is 0 outside)
  arma::mat x = input[variable_in];
  if (variable_in != this->variable)
    return arma::zeros<arma::mat>(arma::size(x));
  arma::mat g = -(x - this->mean)/(this->sd*this->sd);
  g.elem(arma::find((x < this->lower) || (x > this->upper))).fill(arma::datum::nan);
  return g;
}

}
