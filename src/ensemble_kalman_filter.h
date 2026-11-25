#ifndef ENSEMBLEKALMANFILTER_H
#define ENSEMBLEKALMANFILTER_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include <string>

#include "ensemble_kalman.h"
#include "ilike_header.h"
#include "parameters.h"

namespace ilike
{

class EnsembleKalmanOutput;
class ProposalKernel;
class MoveOutput;

class EnsembleKalmanFilter : public EnsembleKalman
{
  
public:
  
  EnsembleKalmanFilter();
  
  EnsembleKalmanFilter(RandomNumberGenerator* rng_in,
                       size_t* seed_in,
                       Data* data_in,
                       size_t lag_in,
                       const std::string &index_name_in,
                       const std::string &time_name_in,
                       const std::string &time_diff_name_in,
                       size_t first_index_in,
                       size_t last_index_in,
                       size_t predictions_per_update_in,
                       double update_time_step_in,
                       double initial_time_in,
                       size_t number_of_ensemble_members_in,
                       EnsembleShifter* shifter_in,
                       std::shared_ptr<Transform> transform_in,
                       IndependentProposalKernel* prior_in,
                       ProposalKernel* transition_model_in,
                       const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in,
                       bool smcfixed_flag_in,
                       bool sequencer_limit_is_fixed_in,
                       bool parallel_in,
                       size_t grain_size_in,
                       const std::string &results_name_in);
  
  virtual ~EnsembleKalmanFilter();
  
  EnsembleKalmanFilter(const EnsembleKalmanFilter &another);
  
  void operator=(const EnsembleKalmanFilter &another);
  LikelihoodEstimator* duplicate() const;
  EnsembleKalman* ensemble_kalman_duplicate() const;

  EnsembleKalmanOutput* specific_ensemble_kalman_initialise();
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters);

  void evaluate(EnsembleKalmanOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(EnsembleKalmanOutput* simulation,
                          const Parameters &conditioned_on_parameters);

  EnsembleKalmanOutput* specific_run();
  EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters);
  
  void point_at_current_schedule(Particle &ensemble_member) const; // (C277)
  // (C278) the update with the current observation: its log-likelihood
  // increment; with commit, the members are shifted (and unpacked), without
  // it (an SMC^2 trial) they are left as they were
  double assimilate(EnsembleKalmanOutput* current_state, bool commit);
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  bool reuses_member_objects() const { return this->predictions_per_update == 1; }
  // (C279) compact ensembles (dev/plan.md, "members stored compactly")
  void move_compact_member(RandomNumberGenerator &rng,
                           const Ensemble* current_ensemble,
                           size_t i,
                           Ensemble* next_ensemble);

  // (dev/plan.md 3.21 (c2''')) SMC^2 with an ensemble Kalman filter: the
  // caller's data tempering takes in the filter's observations one at a time
  // (or a share of one, with anneal). Each output keeps its ensemble, already
  // moved to the next time; a trial is the next observation's likelihood
  // increment from that ensemble (no random numbers); the caller's move stage
  // commits the observation once it is whole (the update, then the move to
  // the next time); a new output is run through the observations committed.
  void set_follows_caller_target(const std::string &fraction_variable, bool anneal);
  bool advances_at_caller_target() const { return this->follows_caller_target; }
  // (C270) parallel outer chunks: a copy runs in each chunk group, drawing
  // from the group's stream (its measurement estimators too); each output owns
  // its own copy of the filter, which is given the main filter's generator
  // afterwards (EnsembleKalmanOutput::repoint_estimator)
  bool supports_parallel_copies() const { return true; }
  void redirect_random_numbers(RandomNumberGenerator* rng_in, size_t* seed_in);
  size_t observations() const { return this->last_index - this->first_index + 1; }
  void advance_at_caller_target_enk(EnsembleKalmanOutput* current_state,
                                    const Parameters &conditioned_on_parameters);
  // (C267) SMC^2's exchange step for the ensemble size (see
  // LikelihoodEstimator::has_exchange)
  void set_exchange(double below, double multiply, size_t maximum);
  bool has_exchange() const { return this->exchange_below >= 0.0; }
  size_t exchange_size(double acceptance_rate, const Parameters &caller_parameters) const;
  size_t get_number_of_inner_particles() const { return this->number_of_ensemble_members; }
  void set_number_of_inner_particles(size_t number_of_inner_particles_in);
  MoveOutput* move_reusing(RandomNumberGenerator &rng,
                           Particle &particle,
                           MoveOutput* spare);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle);

protected:
  
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);

  bool check_termination() const;
  
  //
  void setup_variables();
  
  friend EnsembleKalmanOutput;

  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  // (C265: no smcfixed_flag of its own here: it hid LikelihoodEstimator::smcfixed_flag, which the outputs read)
  
  // point to the one in the hmm_ensemble_factors required to construct this class
  ProposalKernel* transition_model_kernel;
  
  // stored here
  std::unique_ptr<Index> index;

  std::string index_name;
  std::string time_name;
  std::string time_diff_name;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  double initial_time = 0.0; // (each run starts here, C265)
  size_t current_index;
  bool follows_caller_target = false;
  bool compact_conditions() const; // (C279)
  std::string followed_variable;
  bool followed_anneal = false;
  double exchange_below = -1.0; // (C267: no exchange when negative)
  double exchange_multiply = 2.0;
  size_t exchange_max = 0; // (no limit when 0)
  // (the most observations any output has committed; shared by the copies of
  // this filter, since each output owns its own copy of the filter)
  std::shared_ptr<size_t> followed_steps = std::make_shared<size_t>(0);
  void follow_set_time_step(size_t step);
  double follow_increment(EnsembleKalmanOutput* current_state);
  double follow_trial_value(EnsembleKalmanOutput* current_state, const Parameters &conditioned_on_parameters,
                            size_t k, double r);
  void follow_commit(EnsembleKalmanOutput* current_state, const Parameters &conditioned_on_parameters);
  void follow_from_start(EnsembleKalmanOutput* current_state, const Parameters &conditioned_on_parameters);
  void follow_evaluate(EnsembleKalmanOutput* current_state, const Parameters &conditioned_on_parameters);

  // Stored here.
  //EnsembleKalmanFilterOutput* output;
  
  void make_copy(const EnsembleKalmanFilter &another);
  
};
}

#endif
