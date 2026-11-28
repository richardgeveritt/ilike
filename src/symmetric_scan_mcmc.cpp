#include "symmetric_scan_mcmc.h"
#include "deterministic_scan_standard_mcmc_output.h"
#include "particle.h"

namespace ilike
{

SymmetricScanMCMC::SymmetricScanMCMC()
  :DeterministicScanMCMC()
{
}

SymmetricScanMCMC::SymmetricScanMCMC(MCMCTermination* termination_in, const std::vector<MCMC*> &moves_in)
  :DeterministicScanMCMC(termination_in, moves_in)
{
}

SymmetricScanMCMC::~SymmetricScanMCMC()
{
}

SymmetricScanMCMC::SymmetricScanMCMC(const SymmetricScanMCMC &another)
  :DeterministicScanMCMC(another)
{
  // (the same cycle's moves: a copy for another chain is linked again by its cycle)
  this->linked_moves = another.linked_moves;
  this->caller_tuning = another.caller_tuning;
}

void SymmetricScanMCMC::link(const std::vector<const MCMC*> &moves_in)
{
  for (const MCMC* m : moves_in)
    if (m == NULL)
      Rcpp::stop("SymmetricScanMCMC::link: a move of the kernel is missing.");
  this->linked_moves = moves_in;
}

std::vector<const ProposalKernel*> SymmetricScanMCMC::get_proposals() const
{
  if (this->linked_moves.empty())
    return DeterministicScanMCMC::get_proposals();
  std::vector<const ProposalKernel*> proposals;
  for (const MCMC* m : this->linked_moves)
  {
    std::vector<const ProposalKernel*> some = m->get_proposals();
    proposals.insert(proposals.end(), some.begin(), some.end());
  }
  return proposals;
}

// in order, then reversed (each move's reverse), with the caller's tuning
void SymmetricScanMCMC::forward(RandomNumberGenerator &rng, Particle &particle, Particle* scratch) const
{
  if (this->linked_moves.empty())
  {
    if (scratch != NULL)
      DeterministicScanMCMC::move_in_place(rng, particle, *scratch);
    else
      DeterministicScanMCMC::move_in_place(rng, particle);
    return;
  }
  for (const MCMC* m : this->linked_moves)
  {
    if (scratch != NULL)
      m->move_in_place(rng, particle, *scratch);
    else
      m->move_in_place(rng, particle);
  }
}

void SymmetricScanMCMC::backward(RandomNumberGenerator &rng, Particle &particle) const
{
  if (this->linked_moves.empty())
  {
    DeterministicScanMCMC::move_in_place_reversed(rng, particle);
    return;
  }
  for (auto m = this->linked_moves.rbegin(); m != this->linked_moves.rend(); ++m)
    (*m)->move_in_place_reversed(rng, particle);
}

Kernel* SymmetricScanMCMC::duplicate() const
{
  return new SymmetricScanMCMC(*this);
}

MCMC* SymmetricScanMCMC::mcmc_duplicate() const
{
  return new SymmetricScanMCMC(*this);
}

DeterministicScanMCMC* SymmetricScanMCMC::deterministic_scan_mcmc_duplicate() const
{
  return new SymmetricScanMCMC(*this);
}

Particle SymmetricScanMCMC::move(RandomNumberGenerator &rng, const Particle &particle) const
{
  Particle moved = particle;
  this->move_in_place(rng, moved);
  return moved;
}

void SymmetricScanMCMC::move_in_place(RandomNumberGenerator &rng, Particle &particle) const
{
  const ParticleTuning* own = particle.tuning;
  if (this->caller_tuning != NULL)
    particle.tuning = this->caller_tuning;
  this->forward(rng, particle, NULL);
  this->backward(rng, particle);
  particle.tuning = own;
}

void SymmetricScanMCMC::move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch) const
{
  const ParticleTuning* own = particle.tuning;
  if (this->caller_tuning != NULL)
    particle.tuning = this->caller_tuning;
  this->forward(rng, particle, &scratch);
  this->backward(rng, particle);
  particle.tuning = own;
}

void SymmetricScanMCMC::adaptive_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch)
{
  // (the kernel stays fixed: see the class comment)
  this->move_in_place(rng, particle, scratch);
}

void SymmetricScanMCMC::diagnostic_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch,
                                                 bool adapting, double* row)
{
  this->move_in_place(rng, particle, scratch);
}

void SymmetricScanMCMC::move_in_place_reversed(RandomNumberGenerator &rng, Particle &particle) const
{
  this->move_in_place(rng, particle);
}

Particle SymmetricScanMCMC::subsample_move(RandomNumberGenerator &rng, const Particle &particle) const
{
  Rcpp::stop("A symmetric scan of moves (the kernel of an annealed path) is not available with subsampling.");
}

StandardMCMCOutput* SymmetricScanMCMC::initialise_mcmc_output() const
{
  // (it does not adapt: the particles' outputs share it, used read-only)
  return new DeterministicScanStandardMCMCOutput(const_cast<SymmetricScanMCMC*>(this), false);
}

}
