#ifndef UNADJUSTEDMCMC_H
#define UNADJUSTEDMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
class EnsembleKalmanOutput;
class StandardMCMCOutput;
class ProposalKernel;

class UnadjustedMCMC : public MCMC
{
  
public:
  
  UnadjustedMCMC();
  
  // Gaussian random walk.
  UnadjustedMCMC(size_t number_of_iterations_in,
                 const std::vector<Parameters> &initial_points_in,
                 const Parameters &proposal_variances);
  
  UnadjustedMCMC(size_t number_of_iterations_in,
                 ProposalKernel* proposal_in);
  
  UnadjustedMCMC(MCMCTermination* termination_in,
                 ProposalKernel* proposal_in);
  
  virtual ~UnadjustedMCMC();
  
  UnadjustedMCMC(const UnadjustedMCMC &another);
  
  void operator=(const UnadjustedMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  UnadjustedMCMC* unadjusted_mcmc_duplicate() const;
  
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
  ProposalKernel* proposal;
  
  void make_copy(const UnadjustedMCMC &another);
  
  // stored here (change to shared pointer so that memory is not duplicated)
  Index* index;
  
};
}

#endif
