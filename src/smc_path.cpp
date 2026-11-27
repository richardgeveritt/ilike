#include "smc_path.h"
#include "smc_mcmc_move.h"
#include "smc_output.h"
#include "factors.h"
#include "index.h"
#include "mcmc.h"
#include "particle.h"
#include "move_output.h"
#include "independent_proposal_kernel.h"

namespace ilike
{

SMCPath::SMCPath()
{
}

SMCPath::SMCPath(SMCMCMCMove* smc_in,
                 Factors* draw_factors_in,
                 Index* draw_index_in,
                 MCMC* draw_mcmc_in,
                 const std::vector<std::string> &variables_in,
                 const std::vector<std::string> &normaliser_in)
  : smc(smc_in), draw_factors(draw_factors_in), draw_index(draw_index_in), draw_mcmc(draw_mcmc_in),
    variables(variables_in), normaliser(normaliser_in)
{
}

SMCPath::~SMCPath()
{
}

SMCPath::SMCPath(const SMCPath &another)
{
  if (another.smc != NULL)
    this->smc.reset(static_cast<SMCMCMCMove*>(another.smc->smc_duplicate()));
  if (another.draw_factors != NULL)
    this->draw_factors.reset(another.draw_factors->duplicate());
  if (another.draw_index != NULL)
    this->draw_index.reset(another.draw_index->duplicate());
  if (another.draw_mcmc != NULL)
    this->draw_mcmc.reset(another.draw_mcmc->mcmc_duplicate());
  if (another.base != NULL)
    this->base.reset(another.base->independent_proposal_kernel_duplicate());
  this->variables = another.variables;
  this->normaliser = another.normaliser;
  this->kept = another.kept;
}

void SMCPath::set_base(IndependentProposalKernel* base_in)
{
  this->base.reset(base_in);
}

Parameters SMCPath::draw_base(RandomNumberGenerator &rng, const Parameters &given) const
{
  if (this->base == NULL)
    Rcpp::stop("unnormalised(…) { … }: the base measure's simulator was not built.");
  Parameters drawn = this->base->independent_simulate(rng, given);
  Parameters output;
  for (const std::string &v : this->variables)
    output[v] = drawn[v];
  return output;
}

double SMCPath::log_evidence(RandomNumberGenerator &rng, const Parameters &at) const
{
  if (this->smc == NULL)
    Rcpp::stop("normalise = …: this distribution's inner SMC was not built.");
  RandomNumberGenerator local_rng(rng);
  size_t local_seed = size_t(rng());
  local_rng.seed(local_seed);
  this->smc->redirect_random_numbers(&local_rng, &local_seed);
  Parameters conditioned;
  for (const std::string &name : this->normaliser)
    conditioned[name] = at[name];
  std::unique_ptr<SMCOutput> output(this->smc->run(conditioned));
  return output->log_likelihood;
}

AnnealedPath* SMCPath::duplicate() const
{
  return new SMCPath(*this);
}

double SMCPath::single(RandomNumberGenerator &rng, Parameters &state,
                       const Parameters &start, const Parameters &end) const
{
  std::vector<Parameters> states(1, state);
  double log_ratio = this->population(rng, states, start, end);
  for (const std::string &v : this->variables)
    state[v] = states[0][v];
  return log_ratio;
}

double SMCPath::population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                           const Parameters &start, const Parameters &end,
                           double threshold) const
{
  if (this->smc == NULL)
    Rcpp::stop("engine = smc(…): this distribution's inner SMC was not built.");
  // the run's random numbers: a stream of its own, seeded from rng (the SMC
  // seeds its generators from its seed counter)
  RandomNumberGenerator local_rng(rng);
  size_t local_seed = size_t(rng());
  local_rng.seed(local_seed);
  this->smc->redirect_random_numbers(&local_rng, &local_seed);

  Parameters a, b;
  for (const std::string &name : this->normaliser)
  {
    a[name] = start[name];
    b[name] = end[name];
  }
  this->smc->set_path(this->normaliser, a, b);
  std::vector<Parameters> points;
  points.reserve(states.size());
  for (const Parameters &state : states)
  {
    Parameters point;
    for (const std::string &v : this->variables)
      point[v] = state[v];
    points.push_back(point);
  }
  this->smc->set_initial_particles(points, true);
  // (the run given parameters resets the sequencer to the schedule's start)
  std::unique_ptr<SMCOutput> output(this->smc->run(Parameters()));
  double log_ratio = output->log_likelihood;
  // the states at the end of the path (the last generation's points)
  Particles &last = output->back();
  for (size_t i = 0; (i < states.size()) && (i < last.size()); ++i)
    for (const std::string &v : this->variables)
      states[i][v] = last[i]->back().parameters[v];
  return log_ratio;
}

void SMCPath::draw(RandomNumberGenerator &rng, Parameters &state, const Parameters &at, size_t sweeps) const
{
  if (this->draw_mcmc == NULL)
    Rcpp::stop("simulate = sweeps(m): this distribution's MCMC kernel was not built.");
  Parameters point;
  for (const std::string &v : this->variables)
    point[v] = state[v];
  for (const std::string &name : this->normaliser)
    point[name] = at[name];
  this->draw_factors->setup(point);
  Particle particle(point, this->draw_factors.get(), NULL, NULL);
  particle.evaluate_likelihoods(this->draw_index.get());
  for (size_t k = 0; k < sweeps; ++k)
    this->draw_mcmc->move_in_place(rng, particle);
  for (const std::string &v : this->variables)
    state[v] = particle.parameters[v];
}

}
