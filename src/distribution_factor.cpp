#include "distribution_factor.h"

DistributionFactor::DistributionFactor()
  :Factor()
{
}

DistributionFactor::~DistributionFactor()
{
}

DistributionFactor::DistributionFactor(const DistributionFactor &another)
  :Factor(another)
{
  this->make_copy(another);
}

void DistributionFactor::operator=(const DistributionFactor &another)
{
  if(this == &another)
    return;

  Factor::operator=(another);
  this->make_copy(another);
}

void DistributionFactor::make_copy(const DistributionFactor &another)
{
}

double DistributionFactor::evaluate(const Parameters &input) const
{
  return this->distribution_evaluate(input);
}

arma::mat DistributionFactor::evaluate_gradient(const std::string &variable,
                                                const Parameters &input) const
{
  return this->distribution_evaluate_gradient(variable, input);
}
