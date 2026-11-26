#include "unnormalised_gradient_estimator_output.h"
#include "unnormalised_gradient_estimator.h"

namespace ilike
{

UnnormalisedGradientEstimatorOutput::UnnormalisedGradientEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->estimator = NULL;
  this->drawn = false;
  this->gradient = 0.0;
}

UnnormalisedGradientEstimatorOutput::UnnormalisedGradientEstimatorOutput(UnnormalisedGradientEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->drawn = false;
  this->gradient = 0.0;
}

UnnormalisedGradientEstimatorOutput::~UnnormalisedGradientEstimatorOutput()
{
}

UnnormalisedGradientEstimatorOutput::UnnormalisedGradientEstimatorOutput(const UnnormalisedGradientEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void UnnormalisedGradientEstimatorOutput::operator=(const UnnormalisedGradientEstimatorOutput &another)
{
  if (this == &another)
    return;
  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* UnnormalisedGradientEstimatorOutput::duplicate() const
{
  return new UnnormalisedGradientEstimatorOutput(*this);
}

void UnnormalisedGradientEstimatorOutput::make_copy(const UnnormalisedGradientEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->drawn = another.drawn;
  this->gradient = another.gradient;
}

void UnnormalisedGradientEstimatorOutput::simulate()
{
  this->drawn = false;
}

// (new parameters: the fields are drawn at the first gradient asked for)
void UnnormalisedGradientEstimatorOutput::simulate(const Parameters &parameters)
{
  this->drawn = false;
}

// the value is 0 (the normalising constant is handled by the moves)
void UnnormalisedGradientEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
}

void UnnormalisedGradientEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  this->log_likelihood = 0.0;
}

void UnnormalisedGradientEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the gradient of a normalising constant.");
}

void UnnormalisedGradientEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the gradient of a normalising constant.");
}

void UnnormalisedGradientEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the gradient of a normalising constant.");
}

LikelihoodEstimator* UnnormalisedGradientEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat UnnormalisedGradientEstimatorOutput::get_gradient_of_log(const std::string &variable, const Parameters &x)
{
  arma::mat zero(arma::size(x[variable]), arma::fill::zeros);
  if (variable != this->estimator->parameter)
    return zero;
  if (!this->drawn)
  {
    this->gradient = this->estimator->estimate(x);
    this->drawn = true;
  }
  zero.fill(this->gradient);
  return zero;
}

arma::mat UnnormalisedGradientEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable, const Parameters &x)
{
  return this->get_gradient_of_log(variable, x);
}

void UnnormalisedGradientEstimatorOutput::forget_you_were_already_written_to_file()
{
}

void UnnormalisedGradientEstimatorOutput::close_ofstreams()
{
}

void UnnormalisedGradientEstimatorOutput::write_to_file(const std::string &directory_name, const std::string &index)
{
}

}
