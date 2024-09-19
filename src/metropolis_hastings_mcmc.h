#ifndef METROPOLISHASTINGSMCMC_H
#define METROPOLISHASTINGSMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"
#include "proposal_kernel.h"

namespace ilike
{
class StandardMCMCOutput;

class MetropolisHastingsMCMC : public MCMC
{
  
public:
  
  MetropolisHastingsMCMC();
  
  // Gaussian random walk.
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         const std::string &variable_name_in,
                         const arma::mat &proposal_covariance_in);
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         const std::vector<std::string> &variable_names_in,
                         const std::vector<arma::mat> &proposal_covariances_in);
  
  MetropolisHastingsMCMC(size_t number_of_iterations_in,
                         ProposalKernel* proposal_in);
  
  MetropolisHastingsMCMC(MCMCTermination* termination_in,
                         ProposalKernel* proposal_in);
  
  virtual ~MetropolisHastingsMCMC();
  
  MetropolisHastingsMCMC(const MetropolisHastingsMCMC &another);
  
  void operator=(const MetropolisHastingsMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  MetropolisHastingsMCMC* metropolis_hastings_mcmc_duplicate() const;
  
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
  
  void make_copy(const MetropolisHastingsMCMC &another);
  
  // stored here
  Index* index;
  
};
}

#endif
