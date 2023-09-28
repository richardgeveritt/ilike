#ifndef MCMC_H
#define MCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "kernel.h"
#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"

class SMCOutput;
class MCMCTermination;
class MoveOutput;
class MCMCAdaptor;
class StochasticScanMCMC;
class DeterministicScanMCMC;
class Index;
class StandardMCMCOutput;

class MCMC : public Kernel
{

public:

  MCMC();
  MCMC(size_t number_of_iterations_in);
  MCMC(MCMCTermination* termination_in);
  virtual ~MCMC();

  MCMC(const MCMC &another);

  void operator=(const MCMC &another);
  virtual MCMC* mcmc_duplicate() const=0;
  
  MoveOutput* run(RandomNumberGenerator &rng,
                  const Particle &particle) const;
  
  /*
  MoveOutput* run(RandomNumberGenerator &rng,
                  Particle &particle,
                  const Parameters &conditioned_on_parameters);
  */
  
  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            const Particle &particle) const;
  
  /*
  MoveOutput* subsample_run(RandomNumberGenerator &rng,
                            Particle &particle,
                            const Parameters &conditioned_on_parameters);
  */
  
  // could change to also have MoveOutput*
  // change when needed
  // EnsembleMember and particle should merge
  /*
  EnsembleMember run(RandomNumberGenerator &rng,
                     EnsembleMember &particle);
  
  EnsembleMember run(RandomNumberGenerator &rng,
                     EnsembleMember &particle,
                     const Parameters &conditioned_on_parameters);
  
  EnsembleMember subsample_run(RandomNumberGenerator &rng,
                               EnsembleMember &particle,
                               const Parameters &conditioned_on_parameters);
  */
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
  virtual void set_index(Index* index_in)=0;
  virtual void set_index_if_null(Index* index_in)=0;
  
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);
  
  MCMCTermination* get_duplicated_termination() const;
  
  //virtual void set_termination_parameters(h)=0;

protected:
  
  friend StochasticScanMCMC;
  friend DeterministicScanMCMC;
  
  virtual StandardMCMCOutput* initialise_mcmc_output() const=0;
  
  virtual void specific_mcmc_adapt(const Particle &current_particle,
                                   size_t iteration_counter)=0;
  
  // Stored here.
  MCMCTermination* termination;
  
  // stored here
  //Index* index;

  void make_copy(const MCMC &another);

};

#endif
