#include <iterator>
#include "gamma_distribution_factor.h"
#include "distributions.h"

GammaDistributionFactor::GammaDistributionFactor()
  :DistributionFactor()
{
}

GammaDistributionFactor::GammaDistributionFactor(const std::string &variable_in,
                                                 double shape_in,
                                                 double rate_in)
:DistributionFactor()
{
  this->variable = variable_in;
  this->shape = shape_in;
  this->rate = rate_in;
}

GammaDistributionFactor::~GammaDistributionFactor()
{
}

GammaDistributionFactor::GammaDistributionFactor(const GammaDistributionFactor &another)
  :DistributionFactor(another)
{
  this->make_copy(another);
}

void GammaDistributionFactor::operator=(const GammaDistributionFactor &another)
{
  if(this == &another)
    return;

  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor* GammaDistributionFactor::duplicate() const
{
  return( new GammaDistributionFactor(*this));
}

DistributionFactor* GammaDistributionFactor::distribution_factor_duplicate() const
{
  return( new GammaDistributionFactor(*this));
}

void GammaDistributionFactor::make_copy(const GammaDistributionFactor &another)
{
  this->variable = another.variable;
  this->shape = another.shape;
  this->rate = another.rate;
}

double GammaDistributionFactor::distribution_evaluate(const Parameters &input) const
{
  return dgamma(input[this->variable][0],
                this->shape,
                this->rate);
}

arma::mat GammaDistributionFactor::distribution_evaluate_gradient(const std::string &variable,
                                                                     const Parameters &input) const
{
  stop("GammaDistributionFactor::distribution_evaluate_gradient - not written yet.");
}
