#include "unnormalised_normaliser_estimator_output.h"
#include "unnormalised_normaliser_estimator.h"

namespace ilike
{

UnnormalisedNormaliserEstimatorOutput::UnnormalisedNormaliserEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->estimator = NULL;
  this->log_likelihood_smcfixed_part = 0.0;
  this->pending = false;
}

UnnormalisedNormaliserEstimatorOutput::UnnormalisedNormaliserEstimatorOutput(UnnormalisedNormaliserEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->pending = false;
}

UnnormalisedNormaliserEstimatorOutput::~UnnormalisedNormaliserEstimatorOutput()
{
}

UnnormalisedNormaliserEstimatorOutput::UnnormalisedNormaliserEstimatorOutput(const UnnormalisedNormaliserEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void UnnormalisedNormaliserEstimatorOutput::operator=(const UnnormalisedNormaliserEstimatorOutput &another)
{
  if (this == &another)
    return;
  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* UnnormalisedNormaliserEstimatorOutput::duplicate() const
{
  return new UnnormalisedNormaliserEstimatorOutput(*this);
}

void UnnormalisedNormaliserEstimatorOutput::make_copy(const UnnormalisedNormaliserEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->fields = another.fields;
  this->pending = another.pending;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
}

void UnnormalisedNormaliserEstimatorOutput::simulate()
{
  Rcpp::stop("The estimate of a normalising constant needs the parameters.");
}

// the auxiliary variables, M draws of the field at these parameters
void UnnormalisedNormaliserEstimatorOutput::simulate(const Parameters &parameters)
{
  this->pending = true;
  this->draw_if_pending(parameters);
}

// (dev/plan.md 3.21) for a proposed particle: made at the first evaluation
void UnnormalisedNormaliserEstimatorOutput::simulate_when_needed(const Parameters &parameters)
{
  this->fields.clear();
  this->pending = true;
}

void UnnormalisedNormaliserEstimatorOutput::draw_if_pending(const Parameters &parameters)
{
  if (!this->pending)
    return;
  this->fields.resize(this->estimator->points);
  for (size_t m = 0; m < this->estimator->points; ++m)
    this->fields[m] = this->estimator->draw(parameters);
  this->pending = false;
}

void UnnormalisedNormaliserEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  this->draw_if_pending(parameters);
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->estimator->log_estimate(this->fields, parameters);
}

void UnnormalisedNormaliserEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  this->draw_if_pending(parameters);
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->estimator->log_estimate(this->fields, parameters);
}

void UnnormalisedNormaliserEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the estimate of a normalising constant.");
}

void UnnormalisedNormaliserEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the estimate of a normalising constant.");
}

void UnnormalisedNormaliserEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  Rcpp::stop("Subsampling is not available for the estimate of a normalising constant.");
}

LikelihoodEstimator* UnnormalisedNormaliserEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat UnnormalisedNormaliserEstimatorOutput::get_gradient_of_log(const std::string &variable, const Parameters &x)
{
  Rcpp::stop("Gradients are not available for the estimate of a normalising constant (normaliser = sav(...)), so gradient-based "
             "proposals cannot be used for the parameters it depends on.");
}

arma::mat UnnormalisedNormaliserEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable, const Parameters &x)
{
  return this->get_gradient_of_log(variable, x);
}

// (the estimate is not written to the results: its value is part of the
// target, which is)
void UnnormalisedNormaliserEstimatorOutput::forget_you_were_already_written_to_file()
{
}

void UnnormalisedNormaliserEstimatorOutput::close_ofstreams()
{
}

void UnnormalisedNormaliserEstimatorOutput::write_to_file(const std::string &directory_name, const std::string &index)
{
}

}
