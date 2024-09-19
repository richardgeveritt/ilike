#ifndef METROPOLISMCMC_H
#define METROPOLISMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
class SymmetricProposalKernel;
class EnsembleKalmanOutput;
class StandardMCMCOutput;

class MetropolisMCMC : public MCMC
{
  
public:
  
  MetropolisMCMC();
  
  // Gaussian random walk.
  MetropolisMCMC(size_t number_of_iterations_in,
                 const std::vector<Parameters> &initial_points_in,
                 const Parameters &proposal_variances);
  
  MetropolisMCMC(size_t number_of_iterations_in,
                 SymmetricProposalKernel* proposal_in);
  
  MetropolisMCMC(MCMCTermination* termination_in,
                 SymmetricProposalKernel* proposal_in);
  
  virtual ~MetropolisMCMC();
  
  MetropolisMCMC(const MetropolisMCMC &another);
  
  void operator=(const MetropolisMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  MetropolisMCMC* metropolis_mcmc_duplicate() const;
  
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
  SymmetricProposalKernel* proposal;
  
  void make_copy(const MetropolisMCMC &another);
  
  // stored here (change to shared pointer so that memory is not duplicated)
  Index* index;
  
};
}

#endif
