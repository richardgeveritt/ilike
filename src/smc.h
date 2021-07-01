#ifndef SMC_H
#define SMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "smc_output.h"
#include "model_and_algorithm.h"
#include "likelihood_estimator.h"
#include "particles.h"

class SMC : public LikelihoodEstimator
{
public:

  SMC();
  SMC(const ModelAndAlgorithm &model_and_algorithm_in,
      const Data* data_in);
  SMC(const SMC &another);
  virtual ~SMC(void);

  void operator=(const SMC &another);
  virtual SMC* smc_duplicate() const=0;

  SMCOutput* do_smc();

  LikelihoodEstimatorOutput* simulate(const Parameters &parameters);

protected:

  void make_copy(const SMC &another);

  Particles is_step() const;

  //Parameters single_particle_is_step() const;

  virtual void smc_step()=0;

  virtual void weight_update()=0;

  //virtual double single_particle_weight_update() const=0;

  SMCOutput output;
};

#endif
