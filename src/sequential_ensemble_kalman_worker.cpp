#include "sequential_ensemble_kalman_worker.h"
#include "particle_simulator.h"
#include "smc.h"
#include "likelihood_estimator_output.h"
#include "utils.h"
#include "move_output.h"
#include "index.h"
#include "vector_index.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_shifter.h"
#include "factor_variables.h"
#include "single_point_move_output.h"

namespace ilike
{
//Default constructor.
SequentialEnsembleKalmanWorker::SequentialEnsembleKalmanWorker(void)
:EnsembleKalmanWorker()
{
}

SequentialEnsembleKalmanWorker::SequentialEnsembleKalmanWorker(EnsembleKalman* the_enk_in)
:EnsembleKalmanWorker(the_enk_in)
{
  this->log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_ensemble_members());
}

//Copy constructor for the SequentialEnsembleKalmanWorker class.
SequentialEnsembleKalmanWorker::SequentialEnsembleKalmanWorker(const SequentialEnsembleKalmanWorker &another)
:EnsembleKalmanWorker(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialEnsembleKalmanWorker class.
SequentialEnsembleKalmanWorker::~SequentialEnsembleKalmanWorker()
{
}

void SequentialEnsembleKalmanWorker::operator=(const SequentialEnsembleKalmanWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  EnsembleKalmanWorker::operator=(another);
  this->make_copy(another);
}

EnsembleKalmanWorker* SequentialEnsembleKalmanWorker::duplicate() const
{
  return( new SequentialEnsembleKalmanWorker(*this));
}

void SequentialEnsembleKalmanWorker::make_copy(const SequentialEnsembleKalmanWorker &another)
{
  this->log_unnormalised_incremental_weights = another.log_unnormalised_incremental_weights;
}

void SequentialEnsembleKalmanWorker::shift(Ensemble* ensemble,
                                           double inverse_incremental_temperature)
{
  
  this->the_enk->ensemble_shifter->setup(ensemble,
                                         inverse_incremental_temperature);
  
  std::vector<arma::colvec*> measurements = this->the_enk->ensemble_factors->get_measurements();
  bool uses_terms = this->the_enk->ensemble_shifter->uses_shift_terms();
  std::vector<arma::colvec> shift_terms;

  if (ensemble->compact)
  {
    // (C279) a member's shift terms from its rows: the predicted measurement
    // plus the noise times sqrt(inverse incremental temperature), as a
    // member's outputs give them (get_shift)
    std::vector<size_t> active = this->the_enk->ensemble_factors->active_measurement_factor_indices();
    for (size_t i=0; i<ensemble->size(); ++i)
    {
      shift_terms.clear();
      if (uses_terms)
        for (size_t j : active)
          shift_terms.push_back(this->the_enk->ensemble_factors->compact_shift(j,
                                                                             arma::colvec(ensemble->compact_measurement_states[j].row(i).t()),
                                                                             arma::colvec(ensemble->compact_noise[j].row(i).t()),
                                                                             inverse_incremental_temperature));
      this->the_enk->ensemble_shifter->shift(shift_terms,
                                             ensemble->partially_packed_members_col[i],
                                             measurements,
                                             ensemble->kalman_gains,
                                             ensemble->myys,
                                             inverse_incremental_temperature);
    }
    return;
  }
  
  // parallel part (may not be better in parallel)
  for (size_t i=0;
       i<ensemble->members.size();
       ++i)
  {
    if (uses_terms)
      shift_terms = ensemble->members[i]->back().ensemble_factor_variables->get_shifts(inverse_incremental_temperature);
    this->the_enk->ensemble_shifter->shift(shift_terms,
                                           ensemble->partially_packed_members_col[i],
                                           measurements,
                                           ensemble->kalman_gains,
                                           ensemble->myys,
                                           inverse_incremental_temperature);
  }
}

void SequentialEnsembleKalmanWorker::pack(Ensemble* ensemble)
{
  if (ensemble->compact)
  {
    // (C279) the rows are the packed members already: the columns the shift
    // works on, and the active factors' predicted measurements
    size_t N = ensemble->size();
    ensemble->partially_packed_members_col.resize(N);
    for (size_t i = 0; i < N; ++i)
      ensemble->partially_packed_members_col[i] = ensemble->packed_members.row(i).t();
    ensemble->packed_measurement_states.clear();
    for (size_t j : this->the_enk->ensemble_factors->active_measurement_factor_indices())
      ensemble->packed_measurement_states.push_back(ensemble->compact_measurement_states[j]);
    return;
  }
  ensemble->partially_packed_members_row.clear();
  ensemble->partially_packed_members_row.reserve(ensemble->members.size());
  // (the columns too: shift and unpack index them by member, so packing the
  // same ensemble twice must not append a second copy)
  ensemble->partially_packed_members_col.clear();
  ensemble->partially_packed_members_col.reserve(ensemble->members.size());
  
  ensemble->partially_packed_measurement_states.clear();
  ensemble->partially_packed_measurement_states.reserve(ensemble->members.size());

  // can be done in parallel
  for (auto i=ensemble->members.begin();
       i!=ensemble->members.end();
       ++i)
  {
    arma::colvec packed_parameters = (*i)->back().get_colvec(this->the_enk->packing_instructions.states_names);
    ensemble->partially_packed_members_col.push_back(packed_parameters);
    ensemble->partially_packed_members_row.push_back(arma::conv_to<arma::rowvec>::from(packed_parameters));
    ensemble->partially_packed_measurement_states.push_back((*i)->back().ensemble_factor_variables->get_measurement_states_for_covariance());
  }
  
  // must be serial. Each matrix is filled row by row at its final size: it
  // used to grow by join_cols one member at a time, copying O(N^2) numbers for
  // N members (EnKF profile, dev/plan.md C264); the values are the same.
  size_t N = ensemble->partially_packed_members_row.size();
  if (N > 0)
  {
    ensemble->packed_members.set_size(N, ensemble->partially_packed_members_row[0].n_elem);
    for (size_t i = 0; i < N; ++i)
      ensemble->packed_members.row(i) = ensemble->partially_packed_members_row[i];
    size_t number_of_states = ensemble->partially_packed_measurement_states[0].size();
    ensemble->packed_measurement_states.clear();
    ensemble->packed_measurement_states.reserve(number_of_states);
    for (size_t j = 0; j < number_of_states; ++j)
    {
      arma::mat packed(N, ensemble->partially_packed_measurement_states[0][j].n_elem);
      for (size_t i = 0; i < N; ++i)
        packed.row(i) = ensemble->partially_packed_measurement_states[i][j];
      ensemble->packed_measurement_states.push_back(std::move(packed));
    }
  }
  
}

void SequentialEnsembleKalmanWorker::unpack(Ensemble* ensemble)
{
  if (ensemble->compact)
  {
    // (C279) the shifted columns back into the rows
    for (size_t i = 0; i < ensemble->partially_packed_members_col.size(); ++i)
      ensemble->packed_members.row(i) = ensemble->partially_packed_members_col[i].t();
    return;
  }
  // only need to unpack state variables
  
  // can be done in parallel
  for (size_t i=0;
       i<ensemble->members.size();
       ++i)
  {
    for (size_t j=0;
         j<this->the_enk->packing_instructions.states_names.size();
         ++j)
    {
      // (C279) not through operator[], which wrote into a matrix that a copied
      // member (an output copied when an outer SMC resampled) shares: the
      // copy's state moved with this one's update
      ensemble->members[i]->back().parameters.overwrite_variable(this->the_enk->packing_instructions.states_names[j],
                                                                 ensemble->partially_packed_members_col[i](arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                                      this->the_enk->packing_instructions.states_start_and_end[j].second)));
    }
  }
  
}

void SequentialEnsembleKalmanWorker::unpack_with_predicted(Ensemble* ensemble)
{
  if (ensemble->compact)
  {
    // (C281) the analysis into the rows; the predicted members' rows: their
    // states, and, member by member, their predicted measurements and noise
    // simulated from their states alone (as a predicted member object, made
    // from those states, simulates them)
    size_t N = ensemble->partially_packed_members_col.size();
    for (size_t i = 0; i < N; ++i)
      ensemble->packed_members.row(i) = ensemble->partially_packed_members_col[i].t();
    size_t J = ensemble->compact_measurement_states.size();
    ensemble->compact_predicted_states.set_size(N, ensemble->packed_members.n_cols);
    ensemble->compact_predicted_measurement_states.assign(J, arma::mat());
    ensemble->compact_predicted_noise.assign(J, arma::mat());
    arma::colvec measurement_state, noise;
    for (size_t i = 0; i < N; ++i)
    {
      const arma::colvec &predicted = ensemble->partially_packed_predicted_members_col[i];
      ensemble->compact_predicted_states.row(i) = predicted.t();
      Parameters new_parameters;
      for (size_t j=0; j<this->the_enk->packing_instructions.states_names.size(); ++j)
        new_parameters[this->the_enk->packing_instructions.states_names[j]] = predicted(arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                    this->the_enk->packing_instructions.states_start_and_end[j].second));
      for (size_t k = 0; k < J; ++k)
      {
        ensemble->ensemble_factors->compact_simulate(k, new_parameters, measurement_state, noise);
        if (ensemble->compact_predicted_measurement_states[k].n_rows != N)
        {
          ensemble->compact_predicted_measurement_states[k].set_size(N, measurement_state.n_elem);
          ensemble->compact_predicted_noise[k].set_size(N, noise.n_elem);
        }
        ensemble->compact_predicted_measurement_states[k].row(i) = measurement_state.t();
        ensemble->compact_predicted_noise[k].row(i) = noise.t();
      }
    }
    return;
  }
  // only need to unpack state variables
  
  ensemble->predicted_members.clear();
  ensemble->predicted_members.reserve(ensemble->members.size());
  
  // can be done in parallel
  for (size_t i=0;
       i<ensemble->members.size();
       ++i)
  {
    Parameters new_parameters;
    
    for (size_t j=0;
         j<this->the_enk->packing_instructions.states_names.size();
         ++j)
    {
      // (C279) not through operator[], which wrote into a matrix that a copied
      // member (an output copied when an outer SMC resampled) shares: the
      // copy's state moved with this one's update
      ensemble->members[i]->back().parameters.overwrite_variable(this->the_enk->packing_instructions.states_names[j],
                                                                 ensemble->partially_packed_members_col[i](arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                                      this->the_enk->packing_instructions.states_start_and_end[j].second)));

      //ensemble->predicted_members[i]->back().
      
      new_parameters[this->the_enk->packing_instructions.states_names[j]] = ensemble->partially_packed_predicted_members_col[i](arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                                                           this->the_enk->packing_instructions.states_start_and_end[j].second));
    }
    
    ensemble->predicted_members.emplace_back(new SinglePointMoveOutput(new_parameters,ensemble->members[i]->back().ensemble_factor_variables->get_ensemble_factors()));
  }
  
}

void SequentialEnsembleKalmanWorker::weight(Ensemble* ensemble,
                                            const Index* index,
                                            double inverse_incremental_temperature)
{
  ensemble->precompute_gaussian_covariance(inverse_incremental_temperature);
  if (ensemble->compact)
  {
    // (C280) from the members' rows
    for (size_t i = 0; i < ensemble->size(); ++i)
      this->log_unnormalised_incremental_weights[i] = ensemble->ensemble_factors->compact_ensemble_likelihood_ratios(index,
                                                                                                                   ensemble,
                                                                                                                   i,
                                                                                                                   inverse_incremental_temperature,
                                                                                                                   ensemble->inv_sigma_precomps,
                                                                                                                   ensemble->log_det_precomps);
    return;
  }
  for (size_t i = 0; i < ensemble->size(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = (*ensemble)[i]->back().evaluate_ensemble_likelihood_ratios(index,
                                                                                                               inverse_incremental_temperature,
                                                                                                               ensemble->inv_sigma_precomps,
                                                                                                               ensemble->log_det_precomps);
  }
}

void SequentialEnsembleKalmanWorker::subsample_weight(Ensemble* ensemble,
                                                      const Index* index,
                                                      double inverse_incremental_temperature)
{
  if (ensemble->compact)
    Rcpp::stop("SequentialEnsembleKalmanWorker: weights are not available for a compact ensemble yet.");
  ensemble->precompute_gaussian_covariance(inverse_incremental_temperature);
  for (size_t i = 0; i < ensemble->size(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = (*ensemble)[i]->back().subsample_evaluate_ensemble_likelihood_ratios(index,
                                                                                                                         inverse_incremental_temperature,
                                                                                                                         ensemble->inv_sigma_precomps,
                                                                                                                         ensemble->log_det_precomps);
  }
}

arma::colvec SequentialEnsembleKalmanWorker::get_unnormalised_log_incremental_weights() const
{
  return this->log_unnormalised_incremental_weights;
}

// (dev/plan.md 3.22 U9 R3b) correlated estimates: member i's generator at a
// step, seeded from its seed and the step's tag alone (so that its draws do
// not depend on the other members'), with the measurement estimators drawing
// from it while that member is drawn (its perturbed observations), and drawing
// from the estimator's own generator again afterwards
namespace
{
struct MemberStreams
{
  EnsembleKalman* enk;
  EnsembleFactors* factors;
  size_t local_seed = 0;
  MemberStreams(EnsembleKalman* enk_in, EnsembleFactors* factors_in) : enk(enk_in), factors(factors_in) {}
  RandomNumberGenerator member(const RandomNumberGenerator &like, size_t i)
  {
    RandomNumberGenerator g(like);
    g.seed(uint64_t((*this->enk->member_seeds)[i] + this->enk->member_step * 0x9e3779b97f4a7c15ULL));
    return g;
  }
  void draw_from(RandomNumberGenerator &g) { this->factors->redirect_random_numbers(&g, &this->local_seed); }
  ~MemberStreams() { this->factors->redirect_random_numbers(this->enk->get_rng_pointer(), this->enk->get_seed_pointer()); }
};
}

void SequentialEnsembleKalmanWorker::specific_simulate(Ensemble* next_ensemble,
                                                       const Index* index)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  if (this->the_enk->uses_compact_members())
  {
    // (C279) the same draws, in the same order, written into rows
    this->the_enk->start_compact_ensemble(next_ensemble, NULL, NULL);
    for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
      this->the_enk->simulate_compact_member(local_rng, next_ensemble, i, NULL);
    return;
  }
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    Particle* new_particle = next_ensemble->add_ensemble_member();
    this->the_enk->simulate_ensemble_member(local_rng,
                                            new_particle,
                                            this->the_enk->sequencer.schedule_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::specific_simulate(Ensemble* next_ensemble,
                                                       const Index* index,
                                                       const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  if (this->the_enk->member_seeds != NULL)
  {
    // (U9 R3b) each member from its own stream
    this->the_enk->member_step = 0;
    MemberStreams streams(this->the_enk, this->the_enk->ensemble_factors.get());
    const bool compact = this->the_enk->uses_compact_members();
    if (compact)
      this->the_enk->start_compact_ensemble(next_ensemble, NULL, &conditioned_on_parameters);
    for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
    {
      RandomNumberGenerator g = streams.member(local_rng, i);
      streams.draw_from(g);
      if (compact)
        this->the_enk->simulate_compact_member(g, next_ensemble, i, &conditioned_on_parameters);
      else
        this->the_enk->simulate_ensemble_member(g, next_ensemble->add_ensemble_member(),
                                                this->the_enk->sequencer.schedule_parameters, conditioned_on_parameters);
    }
    return;
  }
  
  if (this->the_enk->uses_compact_members())
  {
    this->the_enk->start_compact_ensemble(next_ensemble, NULL, &conditioned_on_parameters);
    for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
      this->the_enk->simulate_compact_member(local_rng, next_ensemble, i, &conditioned_on_parameters);
    return;
  }
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    Particle* new_particle = next_ensemble->add_ensemble_member();
    this->the_enk->simulate_ensemble_member(local_rng,
                                            new_particle,
                                            this->the_enk->sequencer.schedule_parameters,
                                            conditioned_on_parameters);

  }
  
}

void SequentialEnsembleKalmanWorker::specific_move(Ensemble* next_particles,
                                                   Ensemble* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  if (this->the_enk->member_seeds != NULL)
  {
    // (U9 R3b) each member moved with its own stream at this step
    this->the_enk->member_step = this->the_enk->member_step + 1;
    MemberStreams streams(this->the_enk, this->the_enk->ensemble_factors.get());
    if (current_particles->compact)
      this->the_enk->start_compact_ensemble(next_particles, current_particles, NULL);
    for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
    {
      RandomNumberGenerator g = streams.member(local_rng, i);
      streams.draw_from(g);
      if (current_particles->compact)
      {
        this->the_enk->move_compact_member(g, current_particles, i, next_particles);
        continue;
      }
      MoveOutput* spare = next_particles->take_spare();
      if (spare != NULL)
        next_particles->push_back(this->the_enk->move_reusing(g, (*current_particles)[i]->back(), spare));
      else
        next_particles->push_back(this->the_enk->move(g, (*current_particles)[i]->back()));
    }
    return;
  }
  
  if (current_particles->compact)
  {
    // (C279) each member moved from its row into the next ensemble's rows
    this->the_enk->start_compact_ensemble(next_particles, current_particles, NULL);
    for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
      this->the_enk->move_compact_member(local_rng, current_particles, i, next_particles);
    return;
  }
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    // (C264) in a dropped ensemble's member when there is one
    MoveOutput* spare = next_particles->take_spare();
    if (spare != NULL)
      next_particles->push_back(this->the_enk->move_reusing(local_rng,
                                                            (*current_particles)[i]->back(),
                                                            spare));
    else
      next_particles->push_back(this->the_enk->move(local_rng,
                                                    (*current_particles)[i]->back()));
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_specific_move(Ensemble* next_particles,
                                                             Ensemble* current_particles)
{
  if (current_particles->compact)
    Rcpp::stop("SequentialEnsembleKalmanWorker::subsample_specific_move - not available for a compact ensemble.");
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back()));
    
  }
}

// Set up likelihood estimators with all particles.
//for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator l = this->likelihood_estimator_outputs.begin();
//     l!=this->likelihood_estimator_outputs.end();
//     ++l)
//{
//  // Set up each estimator.
//  // Not done yet.
//}

}
