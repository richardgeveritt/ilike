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
  bool is_estimate() const { return true; } // (dev/plan.md 3.21 (c3))

  // (dev/plan.md 3.22 U9 R3b) correlated pseudo-marginal MCMC with this
  // estimate: each member's random numbers from its own seed (its output's
  // member_seeds), and a proposed output keeping the current seeds but k
  // (refresh_members), so that successive estimates are correlated. While a
  // correlated output runs, member_seeds points at its seeds; NULL otherwise
  // (then the members draw as before). The output holds its own copy of the
  // estimator, so the pointer belongs to that copy.
  void set_correlated(bool correlated_in, size_t refresh_members_in)
  {
    this->correlated = correlated_in;
    this->refresh_members = std::max(refresh_members_in, size_t(1));
  }
  bool correlated = false;
  size_t refresh_members = 1;
  const std::vector<uint64_t>* member_seeds = NULL;
  uint64_t member_step = 0; // (the step's tag: 0 at the initial draw, then one more per move)

  
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
  // (C264) as move, building the moved member in a dropped ensemble's member
  // (spare), when reuses_member_objects(); the default deletes spare and moves
  virtual bool reuses_member_objects() const { return false; }
  // (C279) compact ensembles (dev/plan.md, "members stored compactly"):
  // whether this method keeps its members as rows of matrices; the start of a
  // new compact ensemble (initial, or the next one for a move); one member
  // simulated from the proposal or moved, written into its rows. Not used
  // by default (the members stay objects).
  bool uses_compact_members() const;
  void start_compact_ensemble(Ensemble* next_ensemble,
                              const Ensemble* current_ensemble,
                              const Parameters* conditioned_on_parameters);
  void simulate_compact_member(RandomNumberGenerator &rng,
                               Ensemble* next_ensemble,
                               size_t i,
                               const Parameters* conditioned_on_parameters);
  // each method's move of one member (the filter: the transition; EnKI and
  // MFDS: the state kept, its predicted measurement and noise made again)
  virtual void move_compact_member(RandomNumberGenerator &rng,
                                   const Ensemble* current_ensemble,
                                   size_t i,
                                   Ensemble* next_ensemble) { Rcpp::stop("EnsembleKalman::move_compact_member - not available."); }
  // (dev/plan.md 3.21 (c2''')) the caller's move stage, for a filter that
  // follows the caller's data tempering; by default, nothing
  virtual void advance_at_caller_target_enk(EnsembleKalmanOutput* current_state,
                                            const Parameters &conditioned_on_parameters) {}
  virtual MoveOutput* move_reusing(RandomNumberGenerator &rng,
                                   Particle &particle,
                                   MoveOutput* spare);
  
  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     Particle &particle)=0;

protected:
  
  void specific_change_data(Data* new_data);

  // (C279, C280) compact ensembles: a method's own conditions on top of the
  // shared ones (default: not compact); the decision (-1: not made yet; made
  // once the set-up has found the variables); and the two scratch members
  // a compact ensemble's members go through (not copied: a copy makes its own)
  virtual bool compact_conditions() const { return false; }
  // (C281) whether a transform on the members is allowed (MFDS: drawn
  // members are transformed, written ones transformed back), and whether the
  // scratch member is pointed at the current schedule before a move (C277;
  // MFDS's moved members are its predicted members, which carry only their
  // states, so its member objects do not see the schedule there)
  virtual bool compact_allows_transform() const { return false; }
  virtual bool compact_points_at_schedule() const { return true; }
  mutable int compact_decision = -1;
  Particle compact_current;
  Particle compact_next;
  void compact_member_state(const Ensemble* ensemble, size_t i, Particle &member) const;
  void compact_member_rows(Ensemble* ensemble, size_t i, const Parameters &parameters);
  void compact_from_points(Ensemble* next_ensemble, const Parameters* conditioned_on_parameters);
  
  friend EnsembleKalmanOutput;
  friend EnsembleKalmanWorker;
  friend SequentialEnsembleKalmanWorker;
  friend EnsembleSequencer;
  
  EnsembleSequencer sequencer;
  
  // Stored here.
  std::unique_ptr<EnsembleKalmanWorker> the_worker;
  
  EnsembleKalmanOutput* ensemble_kalman_initialise();
  virtual EnsembleKalmanOutput* specific_ensemble_kalman_initialise()=0;

  virtual void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                        const Parameters &conditioned_on_parameters)=0;
  virtual void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                                  const Parameters &conditioned_on_parameters)=0;

  EnsembleKalmanOutput* ensemble_kalman_initialise(const Parameters &parameters);
  virtual EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters)=0;
  virtual void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation)=0;
  virtual void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                        const Parameters &conditioned_on_parameters)=0;
  
  virtual void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                      const Parameters &conditioned_on_parameters)=0;
  virtual void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                        const Parameters &conditioned_on_parameters)=0;

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
  
  void set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in);
  
  // not stored here
  //Parameters* sequencer_parameters;
  
  bool sequencer_limit_is_fixed;
  
  size_t lag;
  size_t number_of_ensemble_members;

  // stored here
  //std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  
  // stored here
  std::unique_ptr<IndependentProposalKernel> proposal;
  
  // stored here
  std::unique_ptr<EnsembleFactors> ensemble_factors;
  
  // stored here
  std::unique_ptr<EnsembleShifter> ensemble_shifter;

  PackingInstructions packing_instructions;
  
  std::vector<std::string> vector_variables;
  std::vector<std::string> any_variables;
  
  std::vector<size_t> vector_variable_sizes;
  
  std::string results_name;
  
  bool proposed_particles_inputted;
  std::vector<Parameters> initial_ensemble; // not needed, unless initial values provided
  
  std::shared_ptr<Transform> transform;
  
  bool initialised;
  
  double reciprocal_schedule_scale;
  
  // Stored here.
  //EnsembleKalmanOutput* output;
  
  std::shared_ptr<HighFive::File> h5_file;
  std::string h5_file_path;
  
  void make_copy(const EnsembleKalman &another);
  
};
}

#endif
