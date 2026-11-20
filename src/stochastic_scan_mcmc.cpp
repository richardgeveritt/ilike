#include "stochastic_scan_mcmc.h"
#include "distributions.h"
#include "stochastic_scan_standard_mcmc_output.h"

namespace ilike
{
// (C226) the moves' diagnostic columns, prefixed with their position, and
// where each move's columns start in a row
static void scan_diagnostic_layout(const std::vector<std::unique_ptr<MCMC>> &moves,
                                   std::vector<std::string> &names,
                                   std::vector<size_t> &offsets,
                                   std::vector<size_t> &sizes)
{
  names.clear();
  offsets.assign(moves.size(), 0);
  sizes.assign(moves.size(), 0);
  for (size_t k = 0; k < moves.size(); ++k)
  {
    std::vector<std::string> own;
    moves[k]->diagnostic_names(own);
    offsets[k] = names.size();
    sizes[k] = own.size();
    for (const std::string &name : own)
      names.push_back("move" + std::to_string(k + 1) + "." + name);
  }
}

StochasticScanMCMC::StochasticScanMCMC()
:MCMC()
{
}

StochasticScanMCMC::StochasticScanMCMC(const std::vector<MCMC*> &moves_in,
                                       const arma::colvec &unnormalised_probabilities_in)
:MCMC()
{
  // takes ownership
  this->moves.clear();
  for (auto i=moves_in.begin(); i!=moves_in.end(); ++i)
    this->moves.emplace_back(*i);
  this->probabilities = unnormalised_probabilities_in/sum(unnormalised_probabilities_in);
}

StochasticScanMCMC::StochasticScanMCMC(MCMCTermination* termination_in,
                                       const std::vector<MCMC*> &moves_in,
                                       const arma::colvec &unnormalised_probabilities_in)
:MCMC(termination_in)
{
  // takes ownership
  this->moves.clear();
  for (auto i=moves_in.begin(); i!=moves_in.end(); ++i)
    this->moves.emplace_back(*i);
  this->probabilities = unnormalised_probabilities_in/sum(unnormalised_probabilities_in);
}

StochasticScanMCMC::~StochasticScanMCMC()
{
}

//Copy constructor for the StochasticScanMCMC class.
StochasticScanMCMC::StochasticScanMCMC(const StochasticScanMCMC &another)
:MCMC(another)
{
  this->make_copy(another);
}

void StochasticScanMCMC::operator=(const StochasticScanMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->moves.clear();
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* StochasticScanMCMC::duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

MCMC* StochasticScanMCMC::mcmc_duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

StochasticScanMCMC* StochasticScanMCMC::stochastic_scan_mcmc_duplicate() const
{
  return( new StochasticScanMCMC(*this));
}

void StochasticScanMCMC::make_copy(const StochasticScanMCMC &another)
{
  this->moves.resize(0);
  this->moves.reserve(another.moves.size());
  for (auto i=another.moves.begin();
       i!=another.moves.end();
       ++i)
  {
    if (*i!=NULL)
      this->moves.emplace_back((*i)->mcmc_duplicate());
    else
      this->moves.emplace_back(nullptr);
  }
  
  this->probabilities = another.probabilities;
}

Particle StochasticScanMCMC::move(RandomNumberGenerator &rng,
                                  const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->move(rng,
                                                           particle);
}

void StochasticScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                       Particle &particle) const
{
  this->moves[rdis(rng, this->probabilities)]->move_in_place(rng, particle);
}

Particle StochasticScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                            const Particle &particle) const
{
  return this->moves[rdis(rng, this->probabilities)]->subsample_move(rng,
                                                                     particle);
}

void StochasticScanMCMC::smc_adapt(SMCOutput* current_state)
{
  // Adapt probabilities?
  
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->smc_adapt(current_state);
  }
}

void StochasticScanMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  // Adapt probabilities?
  
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->ensemble_adapt(current_state);
  }
}

void StochasticScanMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                             size_t iteration_counter)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->specific_mcmc_adapt(current_particle,
                              iteration_counter);
  }
}

void StochasticScanMCMC::set_index(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index(index_in);
  }
}

void StochasticScanMCMC::set_index_if_null(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index_if_null(index_in);
  }
}

void StochasticScanMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}

bool StochasticScanMCMC::reads_cached_factors() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
  {
    if ((*i)->reads_cached_factors())
      return true;
  }
  return false;
}

std::vector<const ProposalKernel*> StochasticScanMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> all_proposals;
  
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    // Every move contributes its proposals (previously only the first move did,
    // so particles had no ProposalStore for the kernels of later moves).
    std::vector<const ProposalKernel*> next_proposals = (*i)->get_proposals();
    all_proposals.insert(all_proposals.end(), next_proposals.begin(), next_proposals.end());
  }
  
  return all_proposals;
}

std::vector<MCMC*> StochasticScanMCMC::get_duplicate_moves() const
{
  std::vector<MCMC*> duplicate_moves;
  duplicate_moves.reserve(this->moves.size());
  
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    duplicate_moves.push_back((*i)->mcmc_duplicate());
  }
  return duplicate_moves;
}

void StochasticScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle,
                                 Particle &scratch) const
{
  this->moves[rdis(rng, this->probabilities)]->move_in_place(rng, particle, scratch);
}

bool StochasticScanMCMC::adapts() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
  {
    if ((*i)->adapts())
      return true;
  }
  return false;
}

StandardMCMCOutput* StochasticScanMCMC::initialise_mcmc_output() const
{
  // A move that adapts gives each chain its own copy of the kernel; otherwise
  // all the particles' outputs share this one (used read-only by the moves).
  if (this->adapts())
    return new StochasticScanStandardMCMCOutput(this->stochastic_scan_mcmc_duplicate());
  return new StochasticScanStandardMCMCOutput(const_cast<StochasticScanMCMC*>(this), false);
}

void StochasticScanMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                                Particle &particle,
                                                Particle &scratch)
{
  // the move chosen adapts its own kernel, if it adapts (dev/plan.md 3.2)
  size_t chosen = rdis(rng, this->probabilities);
  if (this->moves[chosen]->adapts())
    this->moves[chosen]->adaptive_move_in_place(rng, particle, scratch);
  else
    this->moves[chosen]->move_in_place(rng, particle, scratch);
}

void StochasticScanMCMC::append_adaptation_state(std::vector<std::string> *names,
                                 std::vector<double> *values,
                                 bool final,
                                 const std::string &prefix) const
{
  // each move's columns are prefixed with its position in the sweep
  for (size_t k=0; k<this->moves.size(); ++k)
    this->moves[k]->append_adaptation_state(names, values, final,
                                            prefix + "move" + std::to_string(k+1) + ".");
}


void StochasticScanMCMC::diagnostic_names(std::vector<std::string> &names) const
{
  std::vector<size_t> offsets, sizes;
  scan_diagnostic_layout(this->moves, names, offsets, sizes);
}

bool StochasticScanMCMC::warming_up() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
    if ((*i)->warming_up())
      return true;
  return false;
}

void StochasticScanMCMC::diagnostic_move_in_place(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     Particle &scratch,
                                     bool adapting,
                                     double* row)
{
  // the move chosen, as move_in_place / adaptive_move_in_place; the other
  // moves' columns NaN for this iteration
  std::vector<std::string> names;
  std::vector<size_t> offsets, sizes;
  scan_diagnostic_layout(this->moves, names, offsets, sizes);
  for (size_t c = 0; c < names.size(); ++c)
    row[c] = arma::datum::nan;
  size_t chosen = rdis(rng, this->probabilities);
  MCMC* move = this->moves[chosen].get();
  bool adapts = adapting && move->adapts();
  if (sizes[chosen] > 0)
    move->diagnostic_move_in_place(rng, particle, scratch, adapts, row + offsets[chosen]);
  else if (adapts)
    move->adaptive_move_in_place(rng, particle, scratch);
  else
    move->move_in_place(rng, particle, scratch);
}

}
