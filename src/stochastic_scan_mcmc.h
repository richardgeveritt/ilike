#ifndef STOCHASTICSCANMCMC_H
#define STOCHASTICSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

class StochasticScanMCMC : public MCMC
{

public:

  StochasticScanMCMC();

  virtual ~StochasticScanMCMC();

  StochasticScanMCMC(const StochasticScanMCMC &another);

  void operator=(const StochasticScanMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                Particle &particle) const;
  
  Particle move(RandomNumberGenerator &rng,
                Particle &particle,
                const Parameters &conditioned_on_parameters) const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          Particle &particle,
                          const Parameters &conditioned_on_parameters) const;
  
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
  
protected:
  
  void specific_mcmc_adapt(Particle &current_particle,
                           size_t iteration_counter);
  
  // Stored here.
  std::vector<MCMC*> moves;
  
  arma::colvec probabilities;

  void make_copy(const StochasticScanMCMC &another);

};

#endif