#include "annealed_path_kernel.h"
#include "mcmc.h"
#include "particle.h"
#include "factor_variables.h"
#include "vector_factors.h"
#include "likelihood_estimator.h"

namespace ilike
{

AnnealedPathKernel::AnnealedPathKernel(size_t bridges_in, SMCPath* path_in, Index* state_index_in, const std::vector<size_t> &kernel_in)
  : bridges(bridges_in), path(path_in), state_index(state_index_in), kernel(kernel_in)
{
}

AnnealedPathKernel::AnnealedPathKernel(const AnnealedPathKernel &another)
  : bridges(another.bridges), path(another.path), reverse_path(another.reverse_path), kernel(another.kernel)
{
  if (another.state_index != NULL)
    this->state_index.reset(another.state_index->duplicate());
}

AnnealedPathKernel& AnnealedPathKernel::operator=(const AnnealedPathKernel &another)
{
  if (this == &another)
    return *this;
  this->bridges = another.bridges;
  this->path = another.path;
  this->reverse_path = another.reverse_path;
  this->kernel = another.kernel;
  if (another.state_index != NULL)
    this->state_index.reset(another.state_index->duplicate());
  else
    this->state_index.reset();
  return *this;
}

void AnnealedPathKernel::link_kernel(const std::vector<const MCMC*> &cycle, const std::string &function) const
{
  if (!this->path || this->kernel.empty())
    return;
  std::vector<const MCMC*> moves;
  for (size_t k : this->kernel)
  {
    if (k >= cycle.size())
      Rcpp::stop(function + "(…, kernel = …): the kernel's moves are other moves of the algorithm.");
    moves.push_back(cycle[k]);
  }
  this->path->link_kernel(moves);
  if (this->reverse_path)
    this->reverse_path->link_kernel(moves);
}

void AnnealedPathKernel::refuse_estimates(const Particle &current, const Index* index, const std::string &function) const
{
  if (current.factor_variables == NULL)
    return;
  const VectorFactors* factors = dynamic_cast<const VectorFactors*>(current.factor_variables->get_factors());
  if (factors == NULL)
    return;
  for (auto i = index->begin(); i != index->end(); ++i)
  {
    const LikelihoodEstimator* estimator = factors->get_likelihood_estimator(*i);
    if ((estimator != NULL) && estimator->is_estimate())
      Rcpp::stop(function + "(…, bridges = K): factor " + std::to_string(*i + 1) + ", which reads the variables the move changes, "
                 "is an estimate (e.g. a particle filter); annealed moves with estimated factors are not available yet "
                 "(for state-space models they will be AIS with conditional SMC).");
  }
}

}
