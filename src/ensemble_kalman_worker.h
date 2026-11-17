#ifndef ENSEMBLEKALMANWORKER_H
#define ENSEMBLEKALMANWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble.h"
#include "distributions.h"
#include "ilike_header.h"
#include "proposal_kernel.h"
#include "sequencer.h"

namespace ilike
{

class EnsembleKalman;
class ParticleSimulator;
class LikelihoodEstimator;
class EnsembleSequencer;

class EnsembleKalmanWorker
{
public:
  
  EnsembleKalmanWorker();
  
  EnsembleKalmanWorker(EnsembleKalman* the_smc_in);
  
  virtual ~EnsembleKalmanWorker(void);
  
  EnsembleKalmanWorker(const EnsembleKalmanWorker &another);
  void operator=(const EnsembleKalmanWorker &another);
  virtual EnsembleKalmanWorker* duplicate() const=0;
  
  size_t get_number_of_ensemble_members() const;
  RandomNumberGenerator* get_rng();
  size_t get_seed() const;
  void set_seed(size_t seed_in);
  
  void set_enk(EnsembleKalman* the_enk_in);
  
  void simulate(Ensemble* next_ensemble,
                const Index* index);
  void simulate(Ensemble* next_ensemble,
                const Index* index,
                const Parameters &conditioned_on_parameters);
  
  void move(Ensemble* next_particles,
            Ensemble* current_particles);

  void subsample_move(Ensemble* next_particles,
                      Ensemble* current_particles);

  virtual void shift(Ensemble* ensemble,
                     double inverse_incremental_temperature)=0;
  
  virtual void pack(Ensemble* ensemble)=0;
  virtual void unpack(Ensemble* ensemble)=0;
  virtual void unpack_with_predicted(Ensemble* ensemble)=0;
  
  virtual void weight(Ensemble* ensemble,
                      const Index* index,
                      double incremental_temperature)=0;

  virtual void subsample_weight(Ensemble* ensemble,
                                const Index* index,
                                double incremental_temperature)=0;

  virtual arma::colvec get_unnormalised_log_incremental_weights() const=0;
  
protected:
  
  friend EnsembleSequencer;
  
  virtual void specific_simulate(Ensemble* next_ensemble,
                                 const Index* index)=0;
  virtual void specific_simulate(Ensemble* next_ensemble,
                                 const Index* index,
                                 const Parameters &conditioned_on_parameters)=0;
  
  virtual void specific_move(Ensemble* next_particles,
                             Ensemble* current_particles)=0;

  virtual void subsample_specific_move(Ensemble* next_particles,
                                       Ensemble* current_particles)=0;

  void make_copy(const EnsembleKalmanWorker &another);
  
  friend Sequencer;
  EnsembleKalman* the_enk; // not stored here
  
  // stored in model_and_algorithm
  //std::vector<LikelihoodEstimator*> likelihood_estimators;

  // stored here
  // outer vector is which likelihood estimator
  // inner vector is which particle
  //std::vector< std::vector<LikelihoodEstimatorOutput*> > likelihood_estimator_outputs;
  
  // Add in function to do simulation of parts that cannot be done in parallel.
  
};
}

#endif
