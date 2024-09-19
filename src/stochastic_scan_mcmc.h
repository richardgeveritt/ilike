#ifndef STOCHASTICSCANMCMC_H
#define STOCHASTICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

namespace ilike
{
class StochasticScanMCMC : public MCMC
{
  
public:
  
  StochasticScanMCMC();
  
  StochasticScanMCMC(const std::vector<MCMC*> &moves,
                     const arma::colvec &unnormalised_probabilities_in);
  
  StochasticScanMCMC(MCMCTermination* termination_in,
                     const std::vector<MCMC*> &moves_in,
                     const arma::colvec &unnormalised_probabilities_in);
  
  virtual ~StochasticScanMCMC();
  
  StochasticScanMCMC(const StochasticScanMCMC &another);
  
  void operator=(const StochasticScanMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  StochasticScanMCMC* stochastic_scan_mcmc_duplicate() const;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
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
  
  std::vector<MCMC*> get_duplicate_moves() const;
  
  // Stored here.
  std::vector<MCMC*> moves;
  
  arma::colvec probabilities;
  
  void make_copy(const StochasticScanMCMC &another);
  
};
}

#endif
