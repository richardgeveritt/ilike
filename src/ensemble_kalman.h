#ifndef ENSEMBLEKALMAN_H
#define ENSEMBLEKALMAN_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include <memory>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ensemble_sequencer.h"
#include "packing_instructions.h"

class EnsembleKalmanOutput;
class EnsembleKalmanWorker;
class SequentialEnsembleKalmanWorker;
class MeasurementCovarianceEstimator;
class IndependentProposalKernel;
class EnsembleFactors;
class EnsembleShifter;
class MoveOutput;
class EnsembleSequencer;
class Transform;

#include "ensemble.h"

class EnsembleKalman : public LikelihoodEstimator
{

public:

  EnsembleKalman();

  EnsembleKalman(RandomNumberGenerator* rng_in,
                 size_t* seed_in,
                 Data* data_in,
                 size_t number_of_ensemble_members_in,
                 size_t lag_in,
                 EnsembleShifter* shifter_in,
                 std::shared_ptr<Transform> transform_in,
                 bool smcfixed_flag_in,
                 bool sequencer_limit_is_fixed_in,
                 const std::string &results_name_in);

  virtual ~EnsembleKalman();

  EnsembleKalman(const EnsembleKalman &another);

  void operator=(const EnsembleKalman &another);
  virtual EnsembleKalman* ensemble_kalman_duplicate() const=0;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  LikelihoodEstimatorOutput* initialise();
  
  EnsembleKalmanOutput* run();
  EnsembleKalmanOutput* run(const Parameters &conditioned_on_parameters);

  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);
  void setup_variables();
  void setup_variables(const Parameters &parameters);
  
  void set_packing_instructions();
  
  void simulate_ensemble_member(RandomNumberGenerator &rng,
                                Particle* new_particle,
                                const Parameters &sequencer_parameters) const;
  void simulate_ensemble_member(RandomNumberGenerator &rng,
                                Particle* new_particle,
                                const Parameters &sequencer_parameters,
                                const Parameters &conditioned_on_parameters) const;
  
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle)=0;
  
  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     Particle &particle)=0;
  
  /*
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           Particle &particle,
                           const Parameters &conditioned_on_parameters)=0;
  */
  
  //virtual void weight_for_adapting_sequence(Ensemble &current_particles,
  //                                          double incremental_temperature)=0;

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

protected:
  
  void specific_change_data(Data* new_data);
  
  friend EnsembleKalmanOutput;
  friend EnsembleKalmanWorker;
  friend SequentialEnsembleKalmanWorker;
  friend EnsembleSequencer;
  
  EnsembleSequencer sequencer;
  
  // Stored here.
  EnsembleKalmanWorker* the_worker;
  
  EnsembleKalmanOutput* ensemble_kalman_initialise();
  virtual EnsembleKalmanOutput* specific_ensemble_kalman_initialise()=0;
  
  /*
  virtual void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation)=0;
  */
  
  virtual void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                        const Parameters &conditioned_on_parameters)=0;
  virtual void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                                  const Parameters &conditioned_on_parameters)=0;
  
  /*
  virtual void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation)=0;
  */
  
  EnsembleKalmanOutput* ensemble_kalman_initialise(const Parameters &parameters);
  virtual EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters)=0;
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
  
  //void setup_variables_using_candidate_parameters(const Parameters &candidate_parameters);
  
  void simulate_proposal(EnsembleKalmanOutput* simulation,
                         const Index* index);
  void simulate_proposal(EnsembleKalmanOutput* simulation,
                         const Index* index,
                         const Parameters &conditioned_on_parameters);
  
  virtual EnsembleKalmanOutput* specific_run()=0;
  virtual EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  
  void find_measurement_covariances(EnsembleKalmanOutput* simulation);
  
  void set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in);
  
  // not stored here
  //Parameters* sequencer_parameters;
  
  bool sequencer_limit_is_fixed;
  
  size_t lag;
  size_t number_of_ensemble_members;
  
  //bool likelihood_is_evaluated;
  
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
  
  std::vector<std::string> vector_variables;
  std::vector<std::string> any_variables;
  
  std::vector<size_t> vector_variable_sizes;
  
  std::string results_name;
  
  bool proposed_particles_inputted;
  std::vector<Parameters> initial_ensemble; // not needed, unless initial values provided
  
  std::shared_ptr<Transform> transform;
  //TransformPtr inverse_transform;
  
  bool initialised;
  
  double reciprocal_schedule_scale;

  // Stored here.
  //EnsembleKalmanOutput* output;
  
  std::ofstream log_likelihood_file_stream;
  std::ofstream output_lengths_file_stream;
  std::ofstream vector_variables_file_stream;
  std::ofstream vector_variable_sizes_file_stream;
  std::ofstream incremental_log_likelihood_file_stream;
  std::ofstream ess_file_stream;
  std::ofstream schedule_parameters_file_stream;
  std::ofstream vector_points_file_stream;
  std::ofstream any_points_file_stream;
  std::ofstream time_file_stream;

  void make_copy(const EnsembleKalman &another);

};

#endif
