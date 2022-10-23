#ifndef DENSITYLIKELIHOODESTIMATORWORKER_H
#define DENSITYLIKELIHOODESTIMATORWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"
#include "function_pointers.h"

class DensityLikelihoodEstimator;

class DensityLikelihoodEstimatorWorker
{
public:

  DensityLikelihoodEstimatorWorker();
  
  DensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in);

  virtual ~DensityLikelihoodEstimatorWorker(void);

  DensityLikelihoodEstimatorWorker(const DensityLikelihoodEstimatorWorker &another);
  void operator=(const DensityLikelihoodEstimatorWorker &another);
  virtual DensityLikelihoodEstimatorWorker* duplicate() const=0;

  virtual std::vector<Parameters> get_points() const=0;
  size_t get_number_of_points() const;
  RandomNumberGenerator* get_rng();
  size_t get_seed() const;
  void set_seed(size_t seed_in);

  void simulate();
  void simulate(const Parameters &conditioned_on_parameters);
  
  void subsample_simulate(const Parameters &conditioned_on_parameters);

protected:

  virtual void specific_simulate()=0;
  virtual void specific_simulate(const Parameters &conditioned_on_parameters)=0;
  
  virtual void subsample_specific_simulate(const Parameters &conditioned_on_parameters)=0;

  void make_copy(const DensityLikelihoodEstimatorWorker &another);

  DensityLikelihoodEstimator* the_dle; // not stored here

};

#endif
