#ifndef KERNEL_H
#define KERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file kernel.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */



class SMCOutput;
class EnsembleKalmanOutput;
class Index;

class Kernel
{
  
public:
  
  /**
   * @brief Performs the kernel operation.
   */
  Kernel();
  /**
   * @brief Performs the kernel operation.
   *
   * @param number_of_iterations_in The number of iterations.
   */
  Kernel(size_t number_of_iterations_in);
  /**
   * @brief Performs the ~kernel operation.
   */
  virtual ~Kernel();
  
  /**
   * @brief Performs the kernel operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  Kernel(const Kernel &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const Kernel &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  virtual Kernel* duplicate() const=0;
  
  virtual Particle move(RandomNumberGenerator &rng,
                        const Particle &particle) const=0;
  
  /*
   virtual Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual Particle subsample_move(RandomNumberGenerator &rng,
                                  const Particle &particle) const=0;
  
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
  /**
   * @brief Performs the ensemble adapt operation.
   *
   * @param current_state The current state.
   */
  virtual void ensemble_adapt(EnsembleKalmanOutput* current_state)=0;
  
  //virtual void mcmc_adapt(const Particle &current_particle,
  //                        size_t iteration_counter)=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const Kernel &another);
  
};
}

#endif
