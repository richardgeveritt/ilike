#include <algorithm>
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
  if (this->smc != NULL)
  {
    this->smc->set_path_run(true);
    this->pool_idle.push_back(this->smc.get());
  }
}

SMCPath::~SMCPath()
{
}

SMCPath::SMCPath(const SMCPath &another)
{
  // (the original only: a copy's pool starts with it)
  if (another.smc != NULL)
  {
    this->smc.reset(static_cast<SMCMCMCMove*>(another.smc->smc_duplicate()));
    this->pool_idle.push_back(this->smc.get());
  }
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

SMCMCMCMove* SMCPath::acquire() const
{
  std::unique_lock<std::mutex> lock(this->pool_mutex);
  if (!this->pool_idle.empty())
  {
    SMCMCMCMove* idle = this->pool_idle.back();
    this->pool_idle.pop_back();
    return idle;
  }
  // every copy is in use: wait for one, and copy it for this run (the pool
  // grows by one, so a thread waits only while the pool is growing)
  this->pool_freed.wait(lock, [this] { return !this->pool_idle.empty(); });
  SMCMCMCMove* idle = this->pool_idle.back();
  this->pool_idle.pop_back();
  lock.unlock();
  SMCMCMCMove* copy = static_cast<SMCMCMCMove*>(idle->smc_duplicate());
  lock.lock();
  this->pool_copies.emplace_back(copy);
  this->pool_idle.push_back(idle);
  lock.unlock();
  this->pool_freed.notify_one();
  return copy;
}

void SMCPath::release(SMCMCMCMove* smc_in) const
{
  {
    std::lock_guard<std::mutex> lock(this->pool_mutex);
    this->pool_idle.push_back(smc_in);
  }
  this->pool_freed.notify_one();
}

void SMCPath::link_kernel(const std::vector<const MCMC*> &moves)
{
  std::lock_guard<std::mutex> lock(this->pool_mutex);
  if (this->smc != NULL)
    this->smc->link_kernel(moves);
  for (auto &copy : this->pool_copies)
    copy->link_kernel(moves);
}

double SMCPath::log_evidence(RandomNumberGenerator &rng, const Parameters &at) const
{
  if (this->smc == NULL)
    Rcpp::stop("normalise = …: this distribution's inner SMC was not built.");
  Lease lease(this);
  RandomNumberGenerator local_rng(rng);
  size_t local_seed = size_t(rng());
  local_rng.seed(local_seed);
  lease.smc->redirect_random_numbers(&local_rng, &local_seed);
  Parameters conditioned;
  for (const std::string &name : this->normaliser)
    conditioned[name] = at[name];
  std::unique_ptr<SMCOutput> output(lease.smc->run(conditioned));
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

// the path's run from the given points, with the parameters `names` moved
// from their values in start to those in end; the points are made the last
// generation's
double SMCPath::run_from(RandomNumberGenerator &rng, std::vector<Parameters> &points,
                         const std::vector<std::string> &names,
                         const Parameters &start, const Parameters &end,
                         const ParticleTuning* tuning) const
{
  if (this->smc == NULL)
    Rcpp::stop("engine = smc(…): this distribution's inner SMC was not built.");
  Lease lease(this);
  SMCMCMCMove* run_smc = lease.smc;
  run_smc->set_caller_tuning(tuning);
  // the run's random numbers: a stream of its own, seeded from rng (the SMC
  // seeds its generators from its seed counter)
  RandomNumberGenerator local_rng(rng);
  size_t local_seed = size_t(rng());
  local_rng.seed(local_seed);
  run_smc->redirect_random_numbers(&local_rng, &local_seed);

  Parameters a, b;
  for (const std::string &name : names)
  {
    a[name] = start[name];
    b[name] = end[name];
  }
  run_smc->set_path(names, a, b);
  run_smc->set_initial_particles(points, true);
  // (the run given parameters resets the sequencer to the schedule's start)
  std::unique_ptr<SMCOutput> output(run_smc->run(Parameters()));
  double log_ratio = output->log_likelihood;
  run_smc->set_caller_tuning(NULL);
  // the states at the end of the path (the last generation's points)
  Particles &last = output->back();
  for (size_t i = 0; (i < points.size()) && (i < last.size()); ++i)
  {
    const Parameters &final_parameters = last[i]->back().parameters;
    for (auto v = points[i].vector_begin(); v != points[i].vector_end(); ++v)
      v->second.first = std::make_shared<arma::mat>(final_parameters[v->first]);
  }
  return log_ratio;
}

double SMCPath::population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                           const Parameters &start, const Parameters &end,
                           double threshold) const
{
  std::vector<Parameters> points;
  points.reserve(states.size());
  for (const Parameters &state : states)
  {
    Parameters point;
    for (const std::string &v : this->variables)
      point[v] = state[v];
    points.push_back(point);
  }
  double log_ratio = this->run_from(rng, points, this->normaliser, start, end);
  for (size_t i = 0; i < states.size(); ++i)
    for (const std::string &v : this->variables)
      states[i][v] = points[i][v];
  return log_ratio;
}

double SMCPath::single_state(RandomNumberGenerator &rng, Parameters &state,
                             const std::vector<std::string> &names,
                             const Parameters &start, const Parameters &end,
                             const ParticleTuning* tuning) const
{
  // the state's variables other than the path's (and the fixed ones, which
  // the SMC's own parameters supply)
  Parameters point;
  for (auto v = state.vector_begin(); v != state.vector_end(); ++v)
    if (!v->second.second && (std::find(names.begin(), names.end(), v->first) == names.end()))
      point[v->first] = *v->second.first;
  std::vector<Parameters> points(1, point);
  double log_ratio = this->run_from(rng, points, names, start, end, tuning);
  for (auto v = points[0].vector_begin(); v != points[0].vector_end(); ++v)
    state.overwrite_variable(v->first, *v->second.first);
  for (const std::string &name : names)
    state.overwrite_variable(name, end[name]);
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
