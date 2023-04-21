#ifndef KERNEL_H
#define KERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"

class SMCOutput;
class EnsembleKalmanOutput;
class Index;

class Kernel
{

public:

  Kernel();
  Kernel(size_t number_of_iterations_in);
  virtual ~Kernel();

  Kernel(const Kernel &another);

  void operator=(const Kernel &another);
  virtual Kernel* duplicate() const=0;
  
  virtual Particle move(RandomNumberGenerator &rng,
                        Particle &particle) const=0;

  /*
  virtual Particle move(RandomNumberGenerator &rng,
                        Particle &particle,
                        const Parameters &conditioned_on_parameters) const=0;
  */
  
  virtual Particle subsample_move(RandomNumberGenerator &rng,
                                  Particle &particle) const=0;
  
  /*
  virtual Particle subsample_move(RandomNumberGenerator &rng,
                                  Particle &particle,
                                  const Parameters &conditioned_on_parameters) const=0;
  */
  
  /*
  virtual EnsembleMember move(RandomNumberGenerator &rng,
                              const Index* index,
                              EnsembleMember &particle) const=0;
  
  virtual EnsembleMember move(RandomNumberGenerator &rng,
                              const Index* index,
                              EnsembleMember &particle,
                              const Parameters &conditioned_on_parameters) const=0;
  
  virtual EnsembleMember subsample_move(RandomNumberGenerator &rng,
                                        const Index* index,
                                        EnsembleMember &particle,
                                        const Parameters &conditioned_on_parameters) const=0;
  */
  
  virtual void smc_adapt(SMCOutput* current_state)=0;
  virtual void ensemble_adapt(EnsembleKalmanOutput* current_state)=0;
  
  virtual void mcmc_adapt(Particle &current_particle,
                          size_t iteration_counter)=0;
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const Kernel &another);

};

#endif
