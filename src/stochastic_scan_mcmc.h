#ifndef STOCHASTICSCANMCMC_H
#define STOCHASTICSCANMCMC_H

#include <RcppArmadillo.h>
#include <memory>
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

  // the moves' diagnostics (C226), each move's columns prefixed moveK.
  void diagnostic_names(std::vector<std::string> &names) const;
  void diagnostic_move_in_place(RandomNumberGenerator &rng,
                                Particle &particle,
                                Particle &scratch,
                                bool adapting,
                                double* row);
  bool warming_up() const;
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;

  void smc_adapt(SMCOutput* current_state);
  
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool reads_cached_factors() const;
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
protected:
  
  void specific_mcmc_adapt(const Particle &current_particle,
                           size_t iteration_counter);
  
  StandardMCMCOutput* initialise_mcmc_output() const;
  
  std::vector<MCMC*> get_duplicate_moves() const;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<MCMC>> moves;
  
  arma::colvec probabilities;
  
  void make_copy(const StochasticScanMCMC &another);
  
};
}

#endif
