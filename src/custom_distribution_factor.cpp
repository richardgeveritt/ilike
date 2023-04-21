#include <iterator>
#include "custom_distribution_factor.h"

CustomDistributionFactor::CustomDistributionFactor()
  :DistributionFactor()
{
}

CustomDistributionFactor::CustomDistributionFactor(EvaluateLogDistributionPtr distribution_in)
:DistributionFactor()
{
  this->distribution = distribution_in;
  this->distribution_gradient = NULL;
}

CustomDistributionFactor::~CustomDistributionFactor()
{
}

CustomDistributionFactor::CustomDistributionFactor(const CustomDistributionFactor &another)
  :DistributionFactor(another)
{
  this->make_copy(another);
}

void CustomDistributionFactor::operator=(const CustomDistributionFactor &another)
{
  if(this == &another)
    return;

  DistributionFactor::operator=(another);
  this->make_copy(another);
}

Factor* CustomDistributionFactor::duplicate() const
{
  return( new CustomDistributionFactor(*this));
}

DistributionFactor* CustomDistributionFactor::distribution_factor_duplicate() const
{
  return( new CustomDistributionFactor(*this));
}

void CustomDistributionFactor::make_copy(const CustomDistributionFactor &another)
{
  this->distribution = another.distribution;
  this->distribution_gradient = another.distribution_gradient;
}

double CustomDistributionFactor::distribution_evaluate(const Parameters &input) const
{
  return this->distribution(input);
}

arma::mat CustomDistributionFactor::distribution_evaluate_gradient(const std::string &variable,
                                                                   const Parameters &input) const
{
  if (this->distribution_gradient!=NULL)
  {
    return this->distribution_gradient(variable,input);
  }
  else
  {
    Rcpp::stop("CustomDistributionFactor::distribution_evaluate_gradient - gradient not set.");
  }
}
