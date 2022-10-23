#ifndef ENSEMBLEKALMAN_H
#define ENSEMBLEKALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"
//#include "ensemble_member.h"

class EnsembleKalmanOutput;
class EnsembleKalmanWorker;
class SequentialEnsembleKalmanWorker;
class MeasurementCovarianceEstimator;
class IndependentProposalKernel;
class EnsembleFactors;
class EnsembleShifter;
class MoveOutput;
class EnsembleSequencer;

#include "ensemble.h"

class PackingInstructions
{
public:
  std::vector<std::string> states_names;
  //std::vector<std::string> measurements_names;
  std::vector<std::pair<size_t,size_t>> states_start_and_end;
  //std::vector<std::pair<size_t,size_t>> measurements_start_and_end;
};

class EnsembleKalman : public LikelihoodEstimator
{

public:

  EnsembleKalman();

  EnsembleKalman(RandomNumberGenerator* rng_in,
                 size_t* seed_in,
                 Data* data_in,
                 bool smcfixed_flag_in,
                 bool sequencer_limit_is_fixed_in);

  virtual ~EnsembleKalman();

  EnsembleKalman(const EnsembleKalman &another);

  void operator=(const EnsembleKalman &another);
  //LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  LikelihoodEstimatorOutput* initialise();
  virtual EnsembleKalmanOutput* ensemble_kalman_initialise()=0;
  
  EnsembleKalmanOutput* run();
  EnsembleKalmanOutput* run(const Parameters &conditioned_on_parameters);

  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  virtual EnsembleKalmanOutput* ensemble_kalman_initialise(const Parameters &parameters)=0;
  
  Particle simulate_ensemble_member(RandomNumberGenerator &rng) const;
  Particle simulate_ensemble_member(RandomNumberGenerator &rng,
                                    const Parameters &conditioned_on_parameters) const;
  
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle)=0;
  
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle,
                           const Parameters &conditioned_on_parameters)=0;
  
  //virtual void weight_for_adapting_sequence(Ensemble &current_particles,
  //                                          double incremental_temperature)=0;

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

protected:
  
  friend EnsembleKalmanOutput;
  friend EnsembleKalmanWorker;
  friend SequentialEnsembleKalmanWorker;
  friend EnsembleSequencer;
  // Stored here.
  EnsembleKalmanWorker* the_worker;
  
  virtual void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                        const Parameters &conditioned_on_parameters)=0;
  
  virtual void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                                  const Parameters &conditioned_on_parameters)=0;
  
  virtual void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation)=0;
  
  virtual void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters)=0;
  
  virtual void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters)=0;
  virtual void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters)=0;
  
  //virtual void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
  //                                        const Parameters &conditioned_on_parameters)=0;
  //virtual void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
  //                                        const Parameters &conditioned_on_parameters)=0;
  
  virtual void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters)=0;
  virtual void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                                  const Parameters &conditioned_on_parameters)=0;
  
  void simulate_proposal(EnsembleKalmanOutput* simulation,
                         const Index* index);
  void simulate_proposal(EnsembleKalmanOutput* simulation,
                         const Index* index,
                         const Parameters &conditioned_on_parameters);
  
  virtual EnsembleKalmanOutput* specific_run()=0;
  virtual EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  
  void find_measurement_covariances(EnsembleKalmanOutput* simulation);
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  bool sequencer_limit_is_fixed;
  
  size_t lag;
  size_t number_of_ensemble_members;
  
  bool likelihood_is_evaluated;
  
  // stored here
  //std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  
  // stored here
  IndependentProposalKernel* proposal;
  
  // stored here
  EnsembleFactors* ensemble_factors;
  
  // stored here
  EnsembleShifter* ensemble_shifter;
  
  //EnsembleKalmanUpdater* updater;
  //EnsembleKalmanPredictor* predictor;

  PackingInstructions packing_instructions;
  //std::vector<EvaluateLogLikelihoodPtr> numerator_llhds;
  //std::vector<EvaluateLogDistributionPtr> numerator_distributions;
  //std::vector<EvaluateLogLikelihoodPtr> denominator_llhds;
  //std::vector<EvaluateLogDistributionPtr> denominator_distributions;

  // Stored here.
  //EnsembleKalmanOutput* output;

  void make_copy(const EnsembleKalman &another);

};

#endif
