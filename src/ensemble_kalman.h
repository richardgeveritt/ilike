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
#include "ensemble.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{
  /**
   * @file ensemble_kalman.h
   * @brief Defines the EnsembleKalmanOutput class.
   *
   * Stores and manages the output produced by EnsembleKalman. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class EnsembleKalmanOutput
   * @brief The ensemble kalman output class.
   */


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

class EnsembleKalman : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the ensemblekalman operation.
   */
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
  
  /**
   * @brief Performs the ~ensemblekalman operation.
   */
  virtual ~EnsembleKalman();
  
  /**
   * @brief Performs the ensemblekalman operation.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  EnsembleKalman(const EnsembleKalman &another);
  
  /**
   * @brief Assignment operator for EnsembleKalmanOutput.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void operator=(const EnsembleKalman &another);
  /**
   * @brief Creates a deep copy and returns it as a ensemble_kalman pointer.
   *
   * @return The result.
   */
  virtual EnsembleKalman* ensemble_kalman_duplicate() const=0;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* run();
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* run(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  /**
   * @brief Performs the setup variables operation.
   */
  void setup_variables();
  /**
   * @brief Performs the setup variables operation.
   *
   * @param parameters The parameters.
   */
  void setup_variables(const Parameters &parameters);
  
  /**
   * @brief Sets the packing instructions.
   */
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
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend EnsembleKalmanOutput;
  friend EnsembleKalmanWorker;
  friend SequentialEnsembleKalmanWorker;
  friend EnsembleSequencer;
  
  /** @brief The sequencer. */
  EnsembleSequencer sequencer;
  
  // Stored here.
  /** @brief The the worker. */
  EnsembleKalmanWorker* the_worker;
  
  /**
   * @brief Performs the ensemble kalman initialise operation.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* ensemble_kalman_initialise();
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @return The result.
   */
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
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  virtual EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters)=0;
  /**
   * @brief Performs the ensemble kalman simulate operation.
   *
   * @param simulation The simulation.
   */
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
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  virtual EnsembleKalmanOutput* specific_run()=0;
  /**
   * @brief Class-specific implementation for run.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  virtual EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Finds measurement covariances.
   *
   * @param simulation The simulation.
   */
  void find_measurement_covariances(EnsembleKalmanOutput* simulation);
  
  /**
   * @brief Sets the reciprocal schedule scale.
   *
   * @param reciprocal_schedule_scale_in The reciprocal schedule scale.
   */
  void set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in);
  
  // not stored here
  //Parameters* sequencer_parameters;
  
  /** @brief The sequencer limit is fixed. */
  bool sequencer_limit_is_fixed;
  
  /** @brief The lag. */
  size_t lag;
  /** @brief The number of ensemble members. */
  size_t number_of_ensemble_members;
  
  //bool likelihood_is_evaluated;
  
  // stored here
  //std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  
  // stored here
  /** @brief The proposal. */
  IndependentProposalKernel* proposal;
  
  // stored here
  /** @brief The ensemble factors. */
  EnsembleFactors* ensemble_factors;
  
  // stored here
  /** @brief The ensemble shifter. */
  EnsembleShifter* ensemble_shifter;
  
  //EnsembleKalmanUpdater* updater;
  //EnsembleKalmanPredictor* predictor;
  
  /** @brief The packing instructions. */
  PackingInstructions packing_instructions;
  //std::vector<EvaluateLogLikelihoodPtr> numerator_llhds;
  //std::vector<EvaluateLogDistributionPtr> numerator_distributions;
  //std::vector<EvaluateLogLikelihoodPtr> denominator_llhds;
  //std::vector<EvaluateLogDistributionPtr> denominator_distributions;
  
  /** @brief The vector variables. */
  std::vector<std::string> vector_variables;
  /** @brief The any variables. */
  std::vector<std::string> any_variables;
  
  /** @brief The vector variable sizes. */
  std::vector<size_t> vector_variable_sizes;
  
  /** @brief The results name. */
  std::string results_name;
  
  /** @brief The proposed particles inputted. */
  bool proposed_particles_inputted;
  /** @brief The initial ensemble. */
  std::vector<Parameters> initial_ensemble; // not needed, unless initial values provided
  
  /** @brief The transform. */
  std::shared_ptr<Transform> transform;
  //TransformPtr inverse_transform;
  
  /** @brief The initialised. */
  bool initialised;
  
  /** @brief The reciprocal schedule scale. */
  double reciprocal_schedule_scale;
  
  // Stored here.
  //EnsembleKalmanOutput* output;
  
  /** @brief HDF5 output file (kept open for the duration of a run). */
  std::shared_ptr<HighFive::File> h5_file;
  /** @brief Path to the HDF5 output file. */
  std::string h5_file_path;
  
  /**
   * @brief Copies the state of another EnsembleKalmanOutput into this object.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void make_copy(const EnsembleKalman &another);
  
};
}

#endif
