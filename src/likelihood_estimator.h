#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
#include "distributions.h"
#include "data_subsampler.h"
#include "parameters.h"

namespace ilike
{

class LikelihoodEstimatorOutput;
class SMCWorker;
class Factors;

class LikelihoodEstimator
{
  
public:
  
  LikelihoodEstimator();
  LikelihoodEstimator(RandomNumberGenerator* rng_in,
                      size_t* seed_in,
                      Data* data_in,
                      const Parameters &algorithm_parameters_in,
                      bool smcfixed_flag_in);
  virtual ~LikelihoodEstimator();
  
  LikelihoodEstimator(const LikelihoodEstimator &another);
  
  void operator=(const LikelihoodEstimator &another);
  virtual LikelihoodEstimator* duplicate() const=0;
  
  // Initial simulate involves simulating any of the random variables that are needed in the estimator (except for cases where we will update thesee rvs iteratively.
  virtual LikelihoodEstimatorOutput* initialise()=0;
  virtual LikelihoodEstimatorOutput* initialise(const Parameters &parameters)=0;
  
  virtual void setup()=0;
  virtual void setup(const Parameters &parameters)=0;
  
  // To be called if we just want the likelihood, without splitting the estimation into multiple steps.
  //double estimate();
  double estimate(const Parameters &parameters);
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  void change_data_with_raw_pointer(Data* new_data);
  
  Data* get_data() const;
  
  Data* get_current_data() const;
  
  bool get_smcfixed_flag() const;

  // True if this estimator's outputs can be kept in a flat generation
  // (dev/plan.md 2.6a, G2): their smcfixed part uses no random numbers (so it
  // may be evaluated as soon as a point is simulated without changing any
  // random numbers) and they give adaptive summaries. Default: false.
  virtual bool supports_flat_generations() const { return false; }

  // For SMC generations kept flat between MCMC moves (dev/plan.md 2.6a, G3):
  // supports_restore_from_summary: an output can be put back, from its
  // adaptive summary, in the state it had after its last evaluation (the
  // output's simulate uses no random numbers, so new output objects can be
  // made for the purpose); supports_reinitialise: an output can be reset to
  // a new one's state (LikelihoodEstimatorOutput::reinitialise). Default:
  // neither.
  virtual bool supports_restore_from_summary() const { return false; }
  // (dev/plan.md 3.21) its outputs may defer their draws to their first
  // evaluation (LikelihoodEstimatorOutput::simulation_pending)
  virtual bool defers_simulation() const { return false; }
  // (dev/plan.md 3.21) an upper bound on the log value of one term of this
  // estimator (an ABC kernel's maximum), +inf if none is known
  virtual double log_maximum(const Parameters &parameters) const { return arma::datum::inf; }
  // (dev/plan.md 3.21 (c2)) its outputs keep an inner state that follows the
  // caller's SMC target, advanced at the caller's move stage
  // (LikelihoodEstimatorOutput::advance_at_caller_target)
  virtual bool advances_at_caller_target() const { return false; }
  // (dev/plan.md 3.21 (c2''), F2; C267) SMC^2's exchange step: a filter whose
  // number of inner particles (or ensemble members) the caller may increase,
  // pf(..., increase = exchange(...)) or enkf(..., increase = exchange(...))
  virtual bool has_exchange() const { return false; }
  virtual size_t exchange_size(double acceptance_rate, const Parameters &caller_parameters) const { return 0; }
  virtual size_t get_number_of_inner_particles() const { return 0; }
  virtual void set_number_of_inner_particles(size_t number_of_inner_particles_in) {}
  // (dev/plan.md 3.21 (c3)) its value is a Monte Carlo estimate (an inner
  // algorithm, a synthetic likelihood, an estimated normalising constant):
  // moves that need the exact target (slice sampling) refuse it
  virtual bool is_estimate() const { return false; }
  // (dev/plan.md 3.22 U9 R2) its estimates may be negative (signed
  // pseudo-marginal): only MCMC uses it, writing each iteration's sign
  virtual bool may_be_negative() const { return false; }
  // (dev/plan.md 3.21 (c5)) whether this estimator's outputs will be written
  // (an SMC's moves then store their chains' iterations, which are written;
  // otherwise only the final states, which are all that is used)
  virtual void set_outputs_written(bool written) {}
  virtual bool supports_reinitialise() const { return false; }

  // True if this estimator's outputs may be simulated and evaluated for
  // different particles at the same time (dev/plan.md 2.6a, G4): they use no
  // random numbers from the shared generator or seed counter, and change no
  // shared state (the estimator, its data). Default: false (estimators that
  // run an algorithm draw from the shared seed counter).
  virtual bool supports_parallel_chunks() const { return false; }

  // Parallel outer chunks with inner algorithms (dev/plan.md 2.6b, G5):
  // supports_parallel_copies: a copy of this estimator (duplicate()), with
  // its random numbers redirected (redirect_random_numbers), may run in one
  // thread while other copies run in others, and the outputs it makes can be
  // repointed to this estimator afterwards
  // (LikelihoodEstimatorOutput::repoint_estimator). Default: true when no
  // copy is needed at all (supports_parallel_chunks).
  // redirect_random_numbers: draw from this generator and seed counter
  // instead (this estimator and every estimator inside it).
  virtual bool supports_parallel_copies() const { return this->supports_parallel_chunks(); }
  virtual void redirect_random_numbers(RandomNumberGenerator* rng_in,
                                       size_t* seed_in);
  // (C270) the generator and seed counter it draws from (to give an output's
  // own copy of a filter the main filter's, after a parallel chunk group)
  RandomNumberGenerator* get_rng_pointer() const { return this->rng; }
  size_t* get_seed_pointer() const { return this->seed; }
  
  Parameters algorithm_parameters;
  
protected:
  
  virtual void specific_change_data(Data* new_data)=0;
  
  friend SMCWorker;
  
  // Not stored here. Stored in "main'.
  Data* data;
  
  // not stored here
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  size_t* seed;
  
  // Not stored here. Stored in "main'.
  DataSubsampler* subsampler;

  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Factors> factors;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  void make_copy(const LikelihoodEstimator &another);

};
}

#endif
