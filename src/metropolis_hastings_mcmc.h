#ifndef METROPOLISHASTINGSMCMC_H
#define METROPOLISHASTINGSMCMC_H

#include <RcppArmadillo.h>
#include <memory>
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
  // (virtual: an adapting move is copied for each chain through it, and a
  // NoUTurnMCMC must stay one, C210)
  virtual MetropolisHastingsMCMC* metropolis_hastings_mcmc_duplicate() const;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;

  void move_in_place(RandomNumberGenerator &rng,
  
                     Particle &particle) const;

  void move_in_place(RandomNumberGenerator &rng,

                     Particle &particle,

                     Particle &scratch) const;

  void adaptive_move_in_place(RandomNumberGenerator &rng,
                              Particle &particle,
                              Particle &scratch);

  void append_adaptation_state(std::vector<std::string> *names,
                               std::vector<double> *values,
                               bool final,
                               const std::string &prefix) const;

  bool adapts() const;

  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;

  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  void set_cut_index(const Index* cut_index_in);
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
protected:
  
  void specific_mcmc_adapt(const Particle &current_particle,
                           size_t iteration_counter);
  
  StandardMCMCOutput* initialise_mcmc_output() const;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<ProposalKernel> proposal;
  
  void make_copy(const MetropolisHastingsMCMC &another);
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Index> index;
  
};
}

#endif
