#include "adaptive_iterations_mcmc_termination.h"
#include "mcmc.h"

namespace ilike
{

AdaptiveIterationsMCMCTermination::AdaptiveIterationsMCMCTermination()
:MCMCTermination(), iterations(std::make_shared<size_t>(1))
{
}

AdaptiveIterationsMCMCTermination::AdaptiveIterationsMCMCTermination(double c_in, size_t maximum_in)
:MCMCTermination(), iterations(std::make_shared<size_t>(1)), c(c_in), maximum(maximum_in)
{
  if (!(c_in > 0.0) || !(c_in < 1.0))
    Rcpp::stop("ilike::adaptive_iterations: c must be in (0, 1).");
  if (maximum_in < 1)
    Rcpp::stop("ilike::adaptive_iterations: max must be at least 1.");
}

AdaptiveIterationsMCMCTermination::~AdaptiveIterationsMCMCTermination()
{
}

AdaptiveIterationsMCMCTermination::AdaptiveIterationsMCMCTermination(const AdaptiveIterationsMCMCTermination &another)
:MCMCTermination(another)
{
  this->make_copy(another);
}

void AdaptiveIterationsMCMCTermination::operator=(const AdaptiveIterationsMCMCTermination &another)
{
  if (this == &another)
    return;
  MCMCTermination::operator=(another);
  this->make_copy(another);
}

MCMCTermination* AdaptiveIterationsMCMCTermination::duplicate() const
{
  return new AdaptiveIterationsMCMCTermination(*this);
}

void AdaptiveIterationsMCMCTermination::make_copy(const AdaptiveIterationsMCMCTermination &another)
{
  this->counter = another.counter;
  this->iterations = another.iterations; // (shared: the SMC sets it once for every chain)
  this->c = another.c;
  this->maximum = another.maximum;
}

bool AdaptiveIterationsMCMCTermination::terminate()
{
  if (*this->counter >= *this->iterations)
  {
    *this->counter = 0;
    return true;
  }
  return false;
}

void AdaptiveIterationsMCMCTermination::set_parameters(StandardMCMCOutput* mcmc_output)
{
  this->counter = mcmc_output->get_iteration_counter_pointer();
}

size_t AdaptiveIterationsMCMCTermination::iterations_for(double acceptance_rate) const
{
  if (!(acceptance_rate > 0.0))
    return this->maximum;
  if (acceptance_rate >= 1.0)
    return 1;
  double R = std::ceil(std::log(this->c)/std::log(1.0 - acceptance_rate));
  if (!(R >= 1.0))
    return 1;
  if (R > double(this->maximum))
    return this->maximum;
  return size_t(R);
}

}
