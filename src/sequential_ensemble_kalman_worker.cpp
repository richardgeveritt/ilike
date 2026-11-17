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
  
  // parallel part (may not be better in parallel)
  for (size_t i=0;
       i<ensemble->members.size();
       ++i)
  {
    this->the_enk->ensemble_shifter->shift(ensemble->members[i]->back().ensemble_factor_variables.get(),
                                           ensemble->partially_packed_members_col[i],
                                           measurements,
                                           ensemble->kalman_gains,
                                           ensemble->myys,
                                           inverse_incremental_temperature);
  }
}

void SequentialEnsembleKalmanWorker::pack(Ensemble* ensemble)
{
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
  
  // must be serial
  for (size_t i=0;
       i<ensemble->partially_packed_members_row.size();
       ++i)
  {
    if (i==0)
    {
      ensemble->packed_members = ensemble->partially_packed_members_row[i];
      ensemble->packed_measurement_states.clear();
      ensemble->packed_measurement_states.reserve(ensemble->partially_packed_measurement_states[i].size());
      for (auto j=ensemble->partially_packed_measurement_states[i].begin();
           j!=ensemble->partially_packed_measurement_states[i].end();
           ++j)
      {
        ensemble->packed_measurement_states.push_back(*j);
      }
    }
    else
    {
      ensemble->packed_members = join_cols(ensemble->packed_members,ensemble->partially_packed_members_row[i]);
      for (size_t j=0;
           j<ensemble->partially_packed_measurement_states[i].size();
           ++j)
      {
        ensemble->packed_measurement_states[j] = join_cols(ensemble->packed_measurement_states[j],ensemble->partially_packed_measurement_states[i][j]);
      }
    }
  }
  
}

void SequentialEnsembleKalmanWorker::unpack(Ensemble* ensemble)
{
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
      ensemble->members[i]->back().parameters[this->the_enk->packing_instructions.states_names[j]] = ensemble->partially_packed_members_col[i](arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                                                                          this->the_enk->packing_instructions.states_start_and_end[j].second));
    }
  }
  
}

void SequentialEnsembleKalmanWorker::unpack_with_predicted(Ensemble* ensemble)
{
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
      ensemble->members[i]->back().parameters[this->the_enk->packing_instructions.states_names[j]] = ensemble->partially_packed_members_col[i](arma::span(this->the_enk->packing_instructions.states_start_and_end[j].first,
                                                                                                                                                          this->the_enk->packing_instructions.states_start_and_end[j].second));

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

void SequentialEnsembleKalmanWorker::specific_simulate(Ensemble* next_ensemble,
                                                       const Index* index)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
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
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back()));
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_specific_move(Ensemble* next_particles,
                                                             Ensemble* current_particles)
{
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
