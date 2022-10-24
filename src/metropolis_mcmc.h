#ifndef METROPOLISMCMC_H
#define METROPOLISMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc.h"

class SymmetricProposalKernel;
class EnsembleKalmanOutput;

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

  virtual ~MetropolisMCMC();

  MetropolisMCMC(const MetropolisMCMC &another);

  void operator=(const MetropolisMCMC &another);
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
  
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
protected:
  
  void specific_mcmc_adapt(Particle &current_particle,
                           size_t iteration_counter);
  
  // stored here
  SymmetricProposalKernel* proposal;

  void make_copy(const MetropolisMCMC &another);
  
  // stored here
  Index* index;

};

#endif