#include "ensemble_kalman.h"
#include "ensemble_kalman_output.h"
#include "move_output.h"
#include "ensemble_kalman_worker.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "independent_proposal_kernel.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_shifter.h"
#include "ensemble_sequencer.h"
#include "factors.h"
#include "transform.h"

namespace ilike
{
EnsembleKalman::EnsembleKalman()
:LikelihoodEstimator()
{
  this->proposal.reset();
  this->ensemble_factors.reset();
  this->ensemble_shifter.reset();
  this->proposed_particles_inputted = false;
  this->transform = NULL;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  
  this->the_worker.reset();
}

EnsembleKalman::EnsembleKalman(RandomNumberGenerator* rng_in,
                               size_t* seed_in,
                               Data* data_in,
                               size_t number_of_ensemble_members_in,
                               size_t lag_in,
                               EnsembleShifter* shifter_in,
                               std::shared_ptr<Transform> transform_in,
                               bool smcfixed_flag_in,
                               bool sequencer_limit_is_fixed_in,
                               const std::string &results_name_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->number_of_ensemble_members = number_of_ensemble_members_in;
  if (this->number_of_ensemble_members<=1)
    Rcpp::stop("EnsembleKalman::EnsembleKalman - number of ensemble members must be more than 1.");
  
  this->lag = lag_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  
  this->proposal.reset();
  this->ensemble_factors.reset();
  this->ensemble_shifter.reset(shifter_in);
  
  this->results_name = results_name_in;
  
  this->proposed_particles_inputted = false;
  
  this->transform = transform_in;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  
  this->the_worker.reset();
}

EnsembleKalman::~EnsembleKalman()
{
  
  this->proposal.reset();
  
  this->ensemble_factors.reset();
  
  this->ensemble_shifter.reset();
  
  this->the_worker.reset();
  
}

//Copy constructor for the EnsembleKalman class.
EnsembleKalman::EnsembleKalman(const EnsembleKalman &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void EnsembleKalman::operator=(const EnsembleKalman &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal.reset();
  
  this->ensemble_factors.reset();
  
  this->ensemble_shifter.reset();
  
  this->the_worker.reset();

  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

void EnsembleKalman::make_copy(const EnsembleKalman &another)
{
  this->number_of_ensemble_members = another.number_of_ensemble_members;
  this->lag = another.lag;
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  
  this->reciprocal_schedule_scale = another.reciprocal_schedule_scale;

  if (another.proposal!=NULL)
    this->proposal.reset(another.proposal->independent_proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.ensemble_factors!=NULL)
    this->ensemble_factors.reset(another.ensemble_factors->duplicate());
  else
    this->ensemble_factors.reset();
  
  if (another.ensemble_shifter!=NULL)
    this->ensemble_shifter.reset(another.ensemble_shifter->duplicate());
  else
    this->ensemble_shifter.reset();
  
  if (another.the_worker!=NULL)
    this->the_worker.reset(another.the_worker->duplicate());
  else
    this->the_worker.reset();
  this->the_worker->set_enk(this);
  
  this->transform = another.transform;
  
  this->packing_instructions = another.packing_instructions;
  this->vector_variables = another.vector_variables;
  this->any_variables = another.any_variables;
  this->vector_variable_sizes = another.vector_variable_sizes;
  this->results_name = another.results_name;
  
  this->sequencer = another.sequencer;
  
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_ensemble = another.initial_ensemble;
  
  this->initialised = another.initialised;
  this->correlated = another.correlated;
  this->refresh_members = another.refresh_members;
  this->member_seeds = NULL;
  this->crank_nicolson_sd = another.crank_nicolson_sd;
  this->cn_normals = NULL;
}

EnsembleKalmanOutput* EnsembleKalman::run()
{
  this->setup();
  return this->specific_run();
}

EnsembleKalmanOutput* EnsembleKalman::run(const Parameters &conditioned_on_parameters)
{
  this->setup(conditioned_on_parameters);
  return this->specific_run(conditioned_on_parameters);
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise()
{
  return this->ensemble_kalman_initialise();
}

EnsembleKalmanOutput* EnsembleKalman::ensemble_kalman_initialise()
{
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  
  return this->specific_ensemble_kalman_initialise();
}

void EnsembleKalman::setup()
{
  this->setup_variables();
}

void EnsembleKalman::setup(const Parameters &parameters)
{
  this->setup_variables(parameters);
}

void EnsembleKalman::setup_variables()
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
  this->any_variables = dummy_parameters.get_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
  
  if (this->ensemble_factors!=NULL)
    this->ensemble_factors->setup(dummy_parameters);
}

void EnsembleKalman::setup_variables(const Parameters &parameters)
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng,
                                                            parameters);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
  this->any_variables = dummy_parameters.get_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  dummy_parameters.merge_with_fixed(parameters);
  dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
  
  if (this->ensemble_factors!=NULL)
    this->ensemble_factors->setup(dummy_parameters);
}

void EnsembleKalman::set_packing_instructions()
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise(const Parameters &conditioned_on_parameters)
{
  return this->ensemble_kalman_initialise(conditioned_on_parameters);
}

EnsembleKalmanOutput* EnsembleKalman::ensemble_kalman_initialise(const Parameters &conditioned_on_parameters)
{
  
  if (this->initialised==false)
  {
    this->setup(conditioned_on_parameters);
    this->initialised = true;
  }
  return this->specific_ensemble_kalman_initialise(conditioned_on_parameters);
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index)
{
  // Simulate from proposal.
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
  
  if (!this->proposed_particles_inputted)
  {
    
    this->the_worker->simulate(next_ensemble,
                               index);
  }
  else if (this->uses_compact_members())
  {
    this->compact_from_points(next_ensemble, NULL); // (C281)
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble,
                         this->ensemble_factors.get(),
                         this->sequencer.schedule_parameters);
  }
  
  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  //Ensemble* next_ensemble = current_state->add_ensemble();
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_ensemble,
                               index,
                               conditioned_on_parameters);
  }
  else if (this->uses_compact_members())
  {
    this->compact_from_points(next_ensemble, &conditioned_on_parameters); // (C281)
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble,
                         this->ensemble_factors.get(),
                         conditioned_on_parameters,
                         this->sequencer.schedule_parameters);
  }
  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle,
                                              const Parameters &sequencer_parameters) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng),
                        this->ensemble_factors.get(),
                        sequencer_parameters);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng)),
                        this->ensemble_factors.get(),
                        sequencer_parameters);
  }

  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors.get());
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle,
                                              const Parameters &sequencer_parameters,
                                              const Parameters &conditioned_on_parameters) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng,
                                                             conditioned_on_parameters),
                        this->ensemble_factors.get(),
                        conditioned_on_parameters,
                        sequencer_parameters);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng,
                                                                                        conditioned_on_parameters)),
                        this->ensemble_factors.get(),
                        conditioned_on_parameters,
                        sequencer_parameters);
  }

  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors.get(), conditioned_on_parameters);
}

void EnsembleKalman::find_measurement_covariances(EnsembleKalmanOutput* simulation)
{
  simulation->back().find_measurement_covariances();
}

void EnsembleKalman::set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in)
{
  this->reciprocal_schedule_scale = reciprocal_schedule_scale_in;
}

void EnsembleKalman::specific_change_data(Data* new_data)
{
  
}

// (C279, C280) Compact ensembles: the members are rows (dev/plan.md,
// "members stored compactly"). Used when the method's own conditions hold
// (compact_conditions) and every measurement factor is a Gaussian data
// model computed directly (a predicted measurement h(x) and a Gaussian noise
// draw), there is no transform, no initial ensemble is given, the members
// hold only the state variables (no other variable, nothing that is not a
// matrix), and the developer switch ILIKE_RICH_ENSEMBLES is not set (C281:
// an initial ensemble may be given, as MFDS gives one) (it keeps
// member objects, to compare the two paths). The draws are the member
// objects' draws, in the same order: the transition from the move's
// generator, member by member, and each member's measurement noise from the
// estimators' generator, factor by factor.
bool EnsembleKalman::uses_compact_members() const
{
  if (this->compact_decision < 0)
  {
    bool compact = (getenv("ILIKE_RICH_ENSEMBLES") == NULL) &&
      this->compact_conditions() &&
      ((this->transform == NULL) || this->compact_allows_transform()) &&
      (this->ensemble_factors != NULL) &&
      this->ensemble_factors->supports_compact_members() &&
      this->any_variables.empty() && !this->vector_variables.empty() &&
      (this->vector_variables == this->packing_instructions.states_names);
    this->compact_decision = compact ? 1 : 0;
  }
  return this->compact_decision == 1;
}

// a member's state from its row, into a member whose other variables (the
// outer parameters, the schedule) are set already
void EnsembleKalman::compact_member_state(const Ensemble* ensemble, size_t i, Particle &member) const
{
  for (size_t j = 0; j < this->packing_instructions.states_names.size(); ++j)
  {
    const std::pair<size_t, size_t> &span = this->packing_instructions.states_start_and_end[j];
    member.parameters[this->packing_instructions.states_names[j]] =
      arma::colvec(ensemble->packed_members.row(i).subvec(span.first, span.second).t());
  }
}

// a member's rows from its variables: its state (packed as pack does), and
// each measurement factor's predicted measurement and noise (as the member's
// outputs simulate them)
void EnsembleKalman::compact_member_rows(Ensemble* ensemble, size_t i, const Parameters &parameters)
{
  arma::colvec state = parameters.get_colvec(this->packing_instructions.states_names);
  if (ensemble->packed_members.n_cols != state.n_elem)
    ensemble->packed_members.set_size(ensemble->packed_members.n_rows, state.n_elem);
  ensemble->packed_members.row(i) = state.t();
  for (size_t j = 0; j < ensemble->compact_measurement_states.size(); ++j)
  {
    arma::colvec measurement_state, noise;
    this->ensemble_factors->compact_simulate(j, parameters, measurement_state, noise);
    arma::mat &states = ensemble->compact_measurement_states[j];
    arma::mat &noises = ensemble->compact_noise[j];
    if (states.n_cols != measurement_state.n_elem)
      states.set_size(states.n_rows, measurement_state.n_elem);
    if (noises.n_cols != noise.n_elem)
      noises.set_size(noises.n_rows, noise.n_elem);
    states.row(i) = measurement_state.t();
    noises.row(i) = noise.t();
  }
}

// a new compact ensemble: initial (the outer parameters given), or the next
// one for a move (they are the current ensemble's); its members' shared
// variables (the outer parameters and the schedule, as member objects merge
// them in); and the first scratch member pointed at them and at the current
// schedule (C277)
void EnsembleKalman::start_compact_ensemble(Ensemble* next_ensemble,
                                                  const Ensemble* current_ensemble,
                                                  const Parameters* conditioned_on_parameters)
{
  size_t N = this->number_of_ensemble_members;
  size_t J = this->ensemble_factors->number_of_measurement_factors();
  next_ensemble->compact = true;
  next_ensemble->members.clear();
  next_ensemble->spare.clear();
  next_ensemble->packed_members.set_size(N, 0);
  next_ensemble->compact_measurement_states.assign(J, arma::mat(N, 0));
  next_ensemble->compact_noise.assign(J, arma::mat(N, 0));
  if (current_ensemble != NULL)
  {
    next_ensemble->compact_fixed_parameters = current_ensemble->compact_fixed_parameters;
  }
  else
  {
    next_ensemble->compact_fixed_parameters = Parameters();
    if (conditioned_on_parameters != NULL)
      next_ensemble->compact_fixed_parameters.merge_with_fixed(*conditioned_on_parameters,
                                                              this->sequencer.schedule_parameters);
    else
      next_ensemble->compact_fixed_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  }
  this->compact_current = Particle();
  this->compact_current.parameters = next_ensemble->compact_fixed_parameters;
  if (this->compact_points_at_schedule())
    this->compact_current.parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  this->compact_next = Particle();
}

// (C281) a compact ensemble from the given initial points (MFDS), as
// Ensemble::setup makes members from them: each point with the outer
// parameters (if any) and the schedule merged in, its predicted measurement
// and noise simulated, point by point
void EnsembleKalman::compact_from_points(Ensemble* next_ensemble,
                                         const Parameters* conditioned_on_parameters)
{
  this->start_compact_ensemble(next_ensemble, NULL, conditioned_on_parameters);
  size_t N = this->initial_ensemble.size();
  next_ensemble->packed_members.set_size(N, 0);
  for (size_t j = 0; j < next_ensemble->compact_measurement_states.size(); ++j)
  {
    next_ensemble->compact_measurement_states[j].set_size(N, 0);
    next_ensemble->compact_noise[j].set_size(N, 0);
  }
  next_ensemble->unnormalised_log_weights = arma::colvec(N);
  for (size_t i = 0; i < N; ++i)
  {
    Parameters point = this->initial_ensemble[i];
    if (conditioned_on_parameters != NULL)
      point.merge_with_fixed(*conditioned_on_parameters);
    point.merge_with_fixed(this->sequencer.schedule_parameters);
    this->compact_member_rows(next_ensemble, i, point);
  }
}

// a member drawn from the proposal (as simulate_ensemble_member draws it,
// without a transform) and its rows
void EnsembleKalman::simulate_compact_member(RandomNumberGenerator &rng,
                                                   Ensemble* next_ensemble,
                                                   size_t i,
                                                   const Parameters* conditioned_on_parameters)
{
  Parameters drawn = (conditioned_on_parameters != NULL) ? this->proposal->independent_simulate(rng, *conditioned_on_parameters)
                                                         : this->proposal->independent_simulate(rng);
  if (this->transform != NULL)
    drawn = this->transform->transform(drawn); // (C281, MFDS: as simulate_ensemble_member)
  for (auto v = drawn.vector_begin(); v != drawn.vector_end(); ++v)
    this->compact_current.parameters[v->first] = *v->second.first;
  this->compact_member_rows(next_ensemble, i, this->compact_current.parameters);
}

// (C264) defined here, where MoveOutput is complete (its virtual destructor)
MoveOutput* EnsembleKalman::move_reusing(RandomNumberGenerator &rng,
                                         Particle &particle,
                                         MoveOutput* spare)
{
  delete spare;
  return this->move(rng, particle);
}

}
