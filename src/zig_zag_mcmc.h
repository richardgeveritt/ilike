#ifndef ZIGZAGMCMC_H
#define ZIGZAGMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"
#include "proposal_kernel.h"

namespace ilike
{
// Will not be needed later...
class StandardMCMCOutput;

class ZigZagMCMC : public MCMC
{
  
public:
  
  ZigZagMCMC();
  
  // Gaussian random walk.
  ZigZagMCMC(size_t number_of_iterations_in,
             const std::vector<Parameters> &initial_points_in,
             const Parameters &proposal_variances);
  
  ZigZagMCMC(size_t number_of_iterations_in,
             ProposalKernel* proposal_in);
  
  virtual ~ZigZagMCMC();
  
  ZigZagMCMC(const ZigZagMCMC &another);
  
  void operator=(const ZigZagMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  
  // Zig zag can be run to produce a Particle, just as any MCMC can be. This way is can be combined with other MCMC moves.
  // When we call "run" on it, as in the base class, it will use a succession of moves and produce a succession of particles, just as any other MCMC move would.
  // Instead of this, we need to introduce a base class for standard MCMC and another for PDMPs.
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  /*
   EnsembleMember move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle) const;
   
   EnsembleMember move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle,
   const Parameters &conditioned_on_parameters) const;
   
   EnsembleMember subsample_move(RandomNumberGenerator &rng,
   const Index* index,
   EnsembleMember &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  // Overrides that just run the PDMP, rather than calling run in the MCMC base class.
  MoveOutput* run(RandomNumberGenerator &rng,
                  const Particle &particle) const; // override?
  
  /*
   MoveOutput* run(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const; // override?
   */
  
  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            const Particle &particle) const;
  
  /*
   MoveOutput* subsample_run(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  void smc_adapt(SMCOutput* current_state);
  
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
protected:
  
  void specific_mcmc_adapt(const Particle &current_particle,
                           size_t iteration_counter);
  
  StandardMCMCOutput* initialise_mcmc_output() const;
  
  // stored here
  ProposalKernel* proposal;
  
  void make_copy(const ZigZagMCMC &another);
  
};
}

#endif
