#include <iterator>
#include "custom_likelihood_factor.h"

namespace ilike
{
CustomLikelihoodFactor::CustomLikelihoodFactor()
  :LikelihoodFactor()
{
  this->likelihood = NULL;
  this->likelihood_gradient = NULL;
}

CustomLikelihoodFactor::CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                                               Data* data_in)
:LikelihoodFactor(data_in)
{
  this->likelihood = likelihood_in;
  this->likelihood_gradient = NULL;
}

CustomLikelihoodFactor::CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                                               EvaluateGradientLogLikelihoodPtr likelihood_gradient_in,
                                               Data* data_in)
:LikelihoodFactor(data_in)
{
  this->likelihood = likelihood_in;
  this->likelihood_gradient = likelihood_gradient_in;
}

CustomLikelihoodFactor::~CustomLikelihoodFactor()
{
}

CustomLikelihoodFactor::CustomLikelihoodFactor(const CustomLikelihoodFactor &another)
  :LikelihoodFactor(another)
{
  this->make_copy(another);
}

void CustomLikelihoodFactor::operator=(const CustomLikelihoodFactor &another)
{
  if(this == &another)
    return;

  LikelihoodFactor::operator=(another);
  this->make_copy(another);
}

Factor* CustomLikelihoodFactor::duplicate() const
{
  return( new CustomLikelihoodFactor(*this));
}

LikelihoodFactor* CustomLikelihoodFactor::likelihood_factor_duplicate() const
{
  return( new CustomLikelihoodFactor(*this));
}

void CustomLikelihoodFactor::make_copy(const CustomLikelihoodFactor &another)
{
  this->likelihood = another.likelihood;
  this->likelihood_gradient = another.likelihood_gradient;
}

double CustomLikelihoodFactor::likelihood_evaluate(const Parameters &input) const
{
  return this->likelihood(input,*this->data);
}

arma::mat CustomLikelihoodFactor::likelihood_evaluate_gradient(const std::string &variable,
                                                               const Parameters &input) const
{
  if (this->likelihood_gradient!=NULL)
  {
    return this->likelihood_gradient(variable,input,*this->data);
  }
  else
  {
    Rcpp::stop("CustomLikelihoodFactor::likelihood_evaluate_gradient - gradient not set.");
  }
}

void CustomLikelihoodFactor::specific_set_data()
{
  
}
}
