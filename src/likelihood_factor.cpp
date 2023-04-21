#include "likelihood_factor.h"

LikelihoodFactor::LikelihoodFactor()
  :Factor()
{
}

LikelihoodFactor::LikelihoodFactor(Data* data_in)
:Factor()
{
  this->data = data_in;
}

LikelihoodFactor::~LikelihoodFactor()
{
}

LikelihoodFactor::LikelihoodFactor(const LikelihoodFactor &another)
  :Factor(another)
{
  this->make_copy(another);
}

void LikelihoodFactor::operator=(const LikelihoodFactor &another)
{
  if(this == &another)
    return;

  Factor::operator=(another);
  this->make_copy(another);
}

void LikelihoodFactor::make_copy(const LikelihoodFactor &another)
{
  this->data = another.data;
}

double LikelihoodFactor::evaluate(const Parameters &input) const
{
  return this->likelihood_evaluate(input);
}

arma::mat LikelihoodFactor::evaluate_gradient(const std::string &variable,
                                              const Parameters &input) const
{
  return this->likelihood_evaluate_gradient(variable, input);
}

void LikelihoodFactor::set_data(Data* data_in)
{
  this->data = data_in;
  this->specific_set_data();
}
