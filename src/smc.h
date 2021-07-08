#ifndef SMC_H
#define SMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "model_and_algorithm.h"
#include "likelihood_estimator.h"
#include "particles.h"
#include "distributions.h"
#include "parameters.h"

class SMCOutput;
class SMCWorker;

class SMC : public LikelihoodEstimator
{
public:

  SMC();
  SMC(RandomNumberGenerator* rng_in,
      size_t* seed_in,
      const Data* data_in,
      size_t number_of_particles_in,
      size_t lag_in,
      size_t lag_proposed_in);
  SMC(const SMC &another);
  virtual ~SMC(void);

  void operator=(const SMC &another);
  virtual SMC* smc_duplicate() const=0;

  SMCOutput* run(const Parameters &parameters);

  LikelihoodEstimatorOutput* initial_simulate(const Parameters &parameters);

protected:

  virtual void smc_update(SMCOutput* current_state)=0;

  void make_copy(const SMC &another);

  Particles is_step() const;

  friend SMCWorker;
  // Stored here.
  SMCWorker* the_worker;

  //Parameters single_particle_is_step() const;

  //virtual void smc_step()=0;

  //virtual void weight_update()=0;

  //virtual double single_particle_weight_update() const=0;

  // Stored here.
  //SMCOutput* output;

  size_t number_of_particles;

  size_t lag;
  size_t lag_proposed;

};

#endif
