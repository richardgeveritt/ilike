#ifndef DETERMINISTICSCANMCMC_H
#define DETERMINISTICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

class DeterministicScanMCMC : public MCMC
{

public:

  DeterministicScanMCMC();
  
  DeterministicScanMCMC(const std::vector<MCMC*> &moves_in);
  
  DeterministicScanMCMC(MCMCTermination* termination_in,
                        const std::vector<MCMC*> &moves_in);

  virtual ~DeterministicScanMCMC();

  DeterministicScanMCMC(const DeterministicScanMCMC &another);

  void operator=(const DeterministicScanMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
  Particle move(RandomNumberGenerator &rng,
                Particle &particle,
                const Parameters &conditioned_on_parameters) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                Particle &particle,
                const Parameters &conditioned_on_parameters) const;
  */
  
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void set_index(Index* index_in);
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  std::vector<ProposalKernel*> get_proposals() const;
  
protected:
  
  void specific_mcmc_adapt(const Particle &current_particle,
                           size_t iteration_counter);
  
  // Stored here.
  std::vector<MCMC*> moves;
  
  std::vector<size_t> order;

  void make_copy(const DeterministicScanMCMC &another);

};

#endif
