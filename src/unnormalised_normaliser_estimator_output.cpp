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
  this->replicate_log_estimates = another.replicate_log_estimates;
  this->replicate_signs = another.replicate_signs;
  this->pending = another.pending;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->block_seeds = another.block_seeds;
  this->refresh_one_block = another.refresh_one_block;
  this->log_nu = another.log_nu;
}

void UnnormalisedNormaliserEstimatorOutput::simulate()
{
  Rcpp::stop("The estimate of a normalising constant needs the parameters.");
}

// the auxiliary variables, M draws of the field at these parameters
void UnnormalisedNormaliserEstimatorOutput::simulate(const Parameters &parameters)
{
  this->block_seeds.clear();
  this->refresh_one_block = false;
  this->pending = true;
  this->draw_if_pending(parameters);
  // (the draws made now, as before replicates: under data tempering they are
  // made when each replicate enters the target)
  if (this->estimator->tempering_variable.empty())
    this->current_log_estimate(parameters);
}

// (dev/plan.md 3.21) for a proposed particle: made at the first evaluation
void UnnormalisedNormaliserEstimatorOutput::simulate_when_needed(const Parameters &parameters)
{
  this->replicate_log_estimates.clear();
  this->replicate_signs.clear();
  this->block_seeds.clear();
  this->refresh_one_block = false;
  this->pending = true;
}

bool UnnormalisedNormaliserEstimatorOutput::carries_over() const
{
  return this->estimator->block_poisson && this->estimator->block_poisson_correlated;
}

void UnnormalisedNormaliserEstimatorOutput::carry_over_from(const LikelihoodEstimatorOutput &previous,
                                                            const Parameters &previous_parameters,
                                                            const Parameters &new_parameters,
                                                            const std::vector<std::string>* changed)
{
  if (!this->carries_over())
    return;
  const UnnormalisedNormaliserEstimatorOutput* source = dynamic_cast<const UnnormalisedNormaliserEstimatorOutput*>(&previous);
  if ((source == NULL) || (source->block_seeds.size() != this->estimator->block_poisson_blocks))
    return;
  // (the estimate is made at the first evaluation, as for fresh draws)
  this->block_seeds = source->block_seeds;
  this->refresh_one_block = true;
}

// The block-Poisson estimate (dev/plan.md 3.22 U9 R3a; 2210.02734,
// Definition 1 and Algorithm 1). Block l's generator, seeded by seed_l,
// draws chi_l ~ Pois(m) and then chi_l SMC estimates Zhat^(h,l)(theta). With
// Zhat_P their mean, nu ~ Exp(Zhat_P) (nu = E/Zhat_P, E ~ Exp(1)), and
// B^(h,l) = -nu Zhat^(h,l), the estimate of exp(-nu Z(theta)) is
//   Lhat_B = prod_l exp(a/lambda + m) prod_{h <= chi_l} (B^(h,l) - a)/(m lambda),
// unbiased for any a, and positive when a is below every B^(h,l). nu is
// drawn from its proposal whenever u is, so the value is the target of
// (nu, u) over the density they were drawn from:
//   log |Lhat_B| - log q(nu | theta, u) = log |Lhat_B| - log Zhat_P + nu Zhat_P,
// and the Metropolis-Hastings ratio of two such values is Algorithm 1's.
double UnnormalisedNormaliserEstimatorOutput::block_poisson_log_estimate(const Parameters &parameters, double &sign)
{
  const UnnormalisedNormaliserEstimator &e = *this->estimator;
  if (!e.tempering_variable.empty())
    Rcpp::stop("debias = block_poisson(…) under data tempering is not available yet.");
  const size_t L = e.block_poisson_blocks;
  const double m = e.block_poisson_mean, a = e.block_poisson_lower;
  // u: fresh seeds, or the current particle's with one block's replaced
  if (this->block_seeds.size() != L)
  {
    this->block_seeds.resize(L);
    for (size_t l = 0; l < L; ++l)
      this->block_seeds[l] = uint64_t((*e.rng)());
  }
  else if (this->refresh_one_block)
  {
    size_t l = std::min(size_t(runif(*e.rng) * double(L)), L - 1);
    this->block_seeds[l] = uint64_t((*e.rng)());
  }
  this->refresh_one_block = false;
  // the estimates of log Z, block by block
  std::vector<std::vector<double>> log_z(L);
  std::vector<double> all;
  RandomNumberGenerator first_block(*e.rng);
  for (size_t l = 0; l < L; ++l)
  {
    RandomNumberGenerator block(*e.rng);
    block.seed(this->block_seeds[l]);
    int chi = rpois(block, m);
    for (int h = 0; h < chi; ++h)
    {
      log_z[l].push_back(e.plugin_log_z(block, parameters));
      all.push_back(log_z[l].back());
    }
    if (l == 0)
      first_block = block;
  }
  // (no block has an estimate: one more from block 1's generator, for Zhat_P only)
  if (all.empty())
    all.push_back(e.plugin_log_z(first_block, parameters));
  const double log_z_p = ilike::log_mean_exp(all);
  if (!std::isfinite(log_z_p))
  {
    sign = 1.0;
    return -arma::datum::inf;
  }
  // nu ~ Exp(Zhat_P)
  const double E = rexp(*e.rng, 1.0);
  this->log_nu = std::log(E) - log_z_p;
  double log_abs = 0.0;
  sign = 1.0;
  for (size_t l = 0; l < L; ++l)
  {
    log_abs += a / double(L) + m;
    for (double lz : log_z[l])
    {
      double term = (-std::exp(this->log_nu + lz) - a) / (m * double(L));
      if (term < 0.0)
        sign = -sign;
      log_abs += std::log(std::abs(term));
    }
  }
  return log_abs - log_z_p + E;
}

void UnnormalisedNormaliserEstimatorOutput::draw_if_pending(const Parameters &parameters)
{
  if (!this->pending)
    return;
  // (new parameters: no replicate's estimate is drawn yet)
  this->replicate_log_estimates.assign(this->estimator->replicates, arma::datum::nan);
  this->replicate_signs.assign(this->estimator->replicates, 1.0);
  this->pending = false;
}

double UnnormalisedNormaliserEstimatorOutput::current_log_estimate(const Parameters &parameters)
{
  this->draw_if_pending(parameters);
  size_t active = this->estimator->active_replicates(parameters);
  double total = 0.0, product_of_signs = 1.0;
  for (size_t r = 0; r < active; ++r)
  {
    if (std::isnan(this->replicate_log_estimates[r]))
    {
      if (this->estimator->block_poisson)
        this->replicate_log_estimates[r] = this->block_poisson_log_estimate(parameters, this->replicate_signs[r]);
      else
        this->replicate_log_estimates[r] = this->estimator->draw_replicate_log_estimate(parameters, &this->replicate_signs[r]);
    }
    total += this->replicate_log_estimates[r];
    product_of_signs *= this->replicate_signs[r];
  }
  // (U9 R2) the sign of the estimate, the product of the replicates' signs
  this->sign = product_of_signs;
  return total;
}

void UnnormalisedNormaliserEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->current_log_estimate(parameters);
}

void UnnormalisedNormaliserEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->current_log_estimate(parameters);
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
