#include "unnormalised_normaliser_estimator.h"
#include "unnormalised_normaliser_estimator_output.h"

namespace ilike
{

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator()
  :LikelihoodEstimator()
{
  this->potential = NULL;
  this->exact = NULL;
  this->sweep = NULL;
  this->sweeps = 0;
  this->points = 1;
}

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator(RandomNumberGenerator* rng_in,
                                                                 size_t* seed_in,
                                                                 Data* data_in,
                                                                 const std::string &variable_in,
                                                                 UnnormalisedPotentialPtr potential_in,
                                                                 UnnormalisedExactFunctionPtr exact_in,
                                                                 UnnormalisedSweepFunctionPtr sweep_in,
                                                                 size_t sweeps_in,
                                                                 size_t points_in,
                                                                 const std::vector<std::string> &normaliser_in,
                                                                 const std::vector<double> &at_in,
                                                                 bool smcfixed_flag_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->variable = variable_in;
  this->potential = potential_in;
  this->exact = exact_in;
  this->sweep = sweep_in;
  this->sweeps = sweeps_in;
  this->points = std::max(points_in, size_t(1));
  this->normaliser = normaliser_in;
  this->at = at_in;
}

UnnormalisedNormaliserEstimator::~UnnormalisedNormaliserEstimator()
{
}

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator(const UnnormalisedNormaliserEstimator &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void UnnormalisedNormaliserEstimator::operator=(const UnnormalisedNormaliserEstimator &another)
{
  if (this == &another)
    return;
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* UnnormalisedNormaliserEstimator::duplicate() const
{
  return new UnnormalisedNormaliserEstimator(*this);
}

void UnnormalisedNormaliserEstimator::make_copy(const UnnormalisedNormaliserEstimator &another)
{
  this->variable = another.variable;
  this->potential = another.potential;
  this->exact = another.exact;
  this->sweep = another.sweep;
  this->sweeps = another.sweeps;
  this->points = another.points;
  this->normaliser = another.normaliser;
  this->at = another.at;
}

LikelihoodEstimatorOutput* UnnormalisedNormaliserEstimator::initialise()
{
  return new UnnormalisedNormaliserEstimatorOutput(this);
}

LikelihoodEstimatorOutput* UnnormalisedNormaliserEstimator::initialise(const Parameters &parameters)
{
  return this->initialise();
}

void UnnormalisedNormaliserEstimator::setup()
{
}

void UnnormalisedNormaliserEstimator::setup(const Parameters &parameters)
{
}

void UnnormalisedNormaliserEstimator::specific_change_data(Data* new_data)
{
}

arma::mat UnnormalisedNormaliserEstimator::draw(const Parameters &parameters) const
{
  if (this->exact != NULL)
  {
    Parameters drawn = this->exact(*this->rng, parameters, *this->current_data);
    return drawn[this->variable];
  }
  // sweeps (approximate): from the field in the parameters if it is there (a
  // latent field), otherwise from the base measure
  Parameters state = parameters.deep_copy();
  for (size_t k = 0; k < std::max(this->sweeps, size_t(1)); ++k)
    this->sweep(*this->rng, state, *this->current_data);
  return state[this->variable];
}

double UnnormalisedNormaliserEstimator::log_estimate(const std::vector<arma::mat> &fields, const Parameters &parameters) const
{
  Parameters at_parameters = parameters.deep_copy();
  for (size_t k = 0; k < this->normaliser.size(); ++k)
    at_parameters[this->normaliser[k]] = arma::mat(1, 1, arma::fill::value(this->at[k]));
  // log-sum-exp of log gamma(u | theta_hat) - log gamma(u | theta)
  // (a draw that could not be made, e.g. exact draws of an ising field at
  // beta < 0: density 0 at these parameters)
  for (const arma::mat &field : fields)
    if (field.has_nan())
      return -arma::datum::inf;
  std::vector<double> terms(fields.size());
  double top = -arma::datum::inf;
  for (size_t m = 0; m < fields.size(); ++m)
  {
    terms[m] = this->potential(fields[m], at_parameters, *this->current_data) -
               this->potential(fields[m], parameters, *this->current_data);
    top = std::max(top, terms[m]);
  }
  if (!std::isfinite(top))
    return top;
  double sum = 0.0;
  for (double t : terms)
    sum += std::exp(t - top);
  return top + std::log(sum / double(fields.size()));
}

}
