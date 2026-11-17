#ifndef METROPOLISMCMC_H
#define METROPOLISMCMC_H

#include <RcppArmadillo.h>
#include <memory>
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
  std::unique_ptr<SymmetricProposalKernel> proposal;
  
  void make_copy(const MetropolisMCMC &another);
  
  // stored here (change to shared pointer so that memory is not duplicated)
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Index> index;
  
};
}

#endif
