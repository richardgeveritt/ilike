#include "deterministic_scan_mcmc.h"
#include "deterministic_scan_standard_mcmc_output.h"

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

DeterministicScanMCMC::DeterministicScanMCMC()
:MCMC()
{
}

DeterministicScanMCMC::DeterministicScanMCMC(const std::vector<MCMC*> &moves_in)
:MCMC()
{
  // takes ownership
  this->moves.clear();
  for (auto i=moves_in.begin(); i!=moves_in.end(); ++i)
    this->moves.emplace_back(*i);
  this->order.clear();
  this->order.reserve(this->moves.size());
  for (size_t i=0; i<this->moves.size(); ++i)
  {
    this->order.push_back(i);
  }
  link_cycle(this->moves);
}

DeterministicScanMCMC::DeterministicScanMCMC(MCMCTermination* termination_in,
                                             const std::vector<MCMC*> &moves_in)
:MCMC(termination_in)
{
  // takes ownership
  this->moves.clear();
  for (auto i=moves_in.begin(); i!=moves_in.end(); ++i)
    this->moves.emplace_back(*i);
  this->order.clear();
  this->order.reserve(this->moves.size());
  for (size_t i=0; i<this->moves.size(); ++i)
  {
    this->order.push_back(i);
  }
  link_cycle(this->moves);
}

DeterministicScanMCMC::~DeterministicScanMCMC()
{
}

//Copy constructor for the DeterministicScanMCMC class.
DeterministicScanMCMC::DeterministicScanMCMC(const DeterministicScanMCMC &another)
:MCMC(another)
{
  this->make_copy(another);
}

void DeterministicScanMCMC::operator=(const DeterministicScanMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->moves.clear();
  this->order.clear();
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* DeterministicScanMCMC::duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

MCMC* DeterministicScanMCMC::mcmc_duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

DeterministicScanMCMC* DeterministicScanMCMC::deterministic_scan_mcmc_duplicate() const
{
  return( new DeterministicScanMCMC(*this));
}

void DeterministicScanMCMC::make_copy(const DeterministicScanMCMC &another)
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
  
  this->order = another.order;
  link_cycle(this->moves); // (A2: an annealed move's kernel is this copy's moves)
}

Particle DeterministicScanMCMC::move(RandomNumberGenerator &rng,
                                     const Particle &particle) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->move(rng,
                                             current_particle);
  }
  return current_particle;
}

void DeterministicScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                          Particle &particle) const
{
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    this->moves[*i]->move_in_place(rng, particle);
  }
}

Particle DeterministicScanMCMC::subsample_move(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  Particle current_particle = particle;
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    current_particle = this->moves[*i]->subsample_move(rng,
                                                       current_particle);
  }
  return current_particle;
}

void DeterministicScanMCMC::smc_adapt(SMCOutput* current_state)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->smc_adapt(current_state);
  }
}

void DeterministicScanMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->ensemble_adapt(current_state);
  }
}

void DeterministicScanMCMC::specific_mcmc_adapt(const Particle &current_particle,
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

void DeterministicScanMCMC::set_index(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index(index_in);
  }
}

void DeterministicScanMCMC::set_index_if_null(Index* index_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_index_if_null(index_in);
  }
}

void DeterministicScanMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->moves.begin();
       i!=this->moves.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}

bool DeterministicScanMCMC::reads_cached_factors() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
  {
    if ((*i)->reads_cached_factors())
      return true;
  }
  return false;
}

std::vector<const ProposalKernel*> DeterministicScanMCMC::get_proposals() const
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

std::vector<MCMC*> DeterministicScanMCMC::get_duplicate_moves() const
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

void DeterministicScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle,
                                 Particle &scratch) const
{
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    this->moves[*i]->move_in_place(rng, particle, scratch);
  }
}

bool DeterministicScanMCMC::adapts() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
  {
    if ((*i)->adapts())
      return true;
  }
  return false;
}

StandardMCMCOutput* DeterministicScanMCMC::initialise_mcmc_output() const
{
  // A move that adapts gives each chain its own copy of the kernel; otherwise
  // all the particles' outputs share this one (used read-only by the moves).
  if (this->adapts())
    return new DeterministicScanStandardMCMCOutput(this->deterministic_scan_mcmc_duplicate());
  return new DeterministicScanStandardMCMCOutput(const_cast<DeterministicScanMCMC*>(this), false);
}

void DeterministicScanMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                                   Particle &particle,
                                                   Particle &scratch)
{
  // each move that adapts adapts its own kernel (dev/plan.md 3.2)
  for (std::vector<size_t>::const_iterator i=this->order.begin();
       i!=this->order.end();
       ++i)
  {
    if (this->moves[*i]->adapts())
      this->moves[*i]->adaptive_move_in_place(rng, particle, scratch);
    else
      this->moves[*i]->move_in_place(rng, particle, scratch);
  }
}

void DeterministicScanMCMC::append_adaptation_state(std::vector<std::string> *names,
                                 std::vector<double> *values,
                                 bool final,
                                 const std::string &prefix) const
{
  // each move's columns are prefixed with its position in the sweep
  for (size_t k=0; k<this->moves.size(); ++k)
    this->moves[k]->append_adaptation_state(names, values, final,
                                            prefix + "move" + std::to_string(k+1) + ".");
}


void DeterministicScanMCMC::diagnostic_names(std::vector<std::string> &names) const
{
  std::vector<size_t> offsets, sizes;
  scan_diagnostic_layout(this->moves, names, offsets, sizes);
}

bool DeterministicScanMCMC::warming_up() const
{
  for (auto i=this->moves.begin(); i!=this->moves.end(); ++i)
    if ((*i)->warming_up())
      return true;
  return false;
}

void DeterministicScanMCMC::diagnostic_move_in_place(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     Particle &scratch,
                                     bool adapting,
                                     double* row)
{
  // the moves in order, as move_in_place / adaptive_move_in_place
  std::vector<std::string> names;
  std::vector<size_t> offsets, sizes;
  scan_diagnostic_layout(this->moves, names, offsets, sizes);
  for (std::vector<size_t>::const_iterator i=this->order.begin(); i!=this->order.end(); ++i)
  {
    MCMC* move = this->moves[*i].get();
    bool adapts = adapting && move->adapts();
    if (sizes[*i] > 0)
      move->diagnostic_move_in_place(rng, particle, scratch, adapts, row + offsets[*i]);
    else if (adapts)
      move->adaptive_move_in_place(rng, particle, scratch);
    else
      move->move_in_place(rng, particle, scratch);
  }
}


void DeterministicScanMCMC::move_in_place_reversed(RandomNumberGenerator &rng, Particle &particle) const
{
  for (auto i = this->order.rbegin(); i != this->order.rend(); ++i)
    this->moves[*i]->move_in_place_reversed(rng, particle);
}

bool DeterministicScanMCMC::reversible_parts(std::string &why) const
{
  for (auto i = this->order.begin(); i != this->order.end(); ++i)
    if (!this->moves[*i]->reversible_parts(why))
      return false;
  return true;
}

}
