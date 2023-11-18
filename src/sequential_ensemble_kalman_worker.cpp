#include "sequential_ensemble_kalman_worker.h"
#include "particle_simulator.h"
#include "smc.h"
#include "likelihood_estimator_output.h"
#include "utils.h"
#include "move_output.h"
#include "index.h"
#include "vector_single_index.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_shifter.h"
#include "factor_variables.h"
#include "single_point_move_output.h"

//Default constructor.
SequentialEnsembleKalmanWorker::SequentialEnsembleKalmanWorker(void)
  :EnsembleKalmanWorker()
{
}

SequentialEnsembleKalmanWorker::SequentialEnsembleKalmanWorker(EnsembleKalman* the_enk_in)
  :EnsembleKalmanWorker(the_enk_in)
{
  //this->particles = NULL;
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
  
  //this->particles.clear();

  EnsembleKalmanWorker::operator=(another);
  this->make_copy(another);
}

EnsembleKalmanWorker* SequentialEnsembleKalmanWorker::duplicate() const
{
  return( new SequentialEnsembleKalmanWorker(*this));
}

/*
EnsembleKalman SequentialEnsembleKalmanWorker::simulated_particles() const
{
  return this->particles;
}

EnsembleKalman& SequentialEnsembleKalmanWorker::simulated_particles()
{
  return this->particles;
}
*/

void SequentialEnsembleKalmanWorker::make_copy(const SequentialEnsembleKalmanWorker &another)
{
  //this->particles = another.particles;
  this->log_unnormalised_incremental_weights = another.log_unnormalised_incremental_weights;
}

void SequentialEnsembleKalmanWorker::shift(Ensemble* ensemble,
                                           double inverse_incremental_temperature)
{
  //double inverse_incremental_temperature = ensemble->get_inverse_incremental_temperature();
  
  this->the_enk->ensemble_shifter->setup(ensemble,
                                         inverse_incremental_temperature);
  
  std::vector<arma::colvec*> measurements = this->the_enk->ensemble_factors->get_measurements();
  
  // parallel part (may not be better in parallel)
  for (size_t i=0;
       i<ensemble->members.size();
       ++i)
  {
    this->the_enk->ensemble_shifter->shift(ensemble->members[i]->back().ensemble_factor_variables,
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
  //ensemble->partially_packed_members_col.clear();
  //ensemble->partially_packed_members_col.reserve(ensemble->members.size());
  
  ensemble->partially_packed_measurement_states.clear();
  ensemble->partially_packed_measurement_states.reserve(ensemble->members.size());
  
  //ensemble->partially_packed_measurement_random_shifts.clear();
  //ensemble->partially_packed_measurement_random_shifts.reserve(ensemble->members.size());
  
  // can be done in parallel
  for (auto i=ensemble->members.begin();
       i!=ensemble->members.end();
       ++i)
  {
    arma::colvec packed_parameters = (*i)->back().get_colvec(this->the_enk->packing_instructions.states_names);
    ensemble->partially_packed_members_col.push_back(packed_parameters);
    ensemble->partially_packed_members_row.push_back(arma::conv_to<arma::rowvec>::from(packed_parameters));
    ensemble->partially_packed_measurement_states.push_back((*i)->back().ensemble_factor_variables->get_measurement_states_for_covariance());
    //ensemble->partially_packed_measurement_random_shifts.push_back(i->ensemble_factor_variables->get_measurement_random_shifts());
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
  
  //ensemble->measurements = ensemble->ensemble_factors->get_measurements();
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
    
    ensemble->predicted_members.push_back(new SinglePointMoveOutput(new_parameters,ensemble->members[i]->back().ensemble_factor_variables->get_ensemble_factors()));
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
                                                                                                               inverse_incremental_temperature);
  }
}

/*
void SequentialEnsembleKalmanWorker::weight(Ensemble* ensemble,
                                            const Index* index,
                                            double inverse_incremental_temperature,
                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < ensemble->size(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = (*ensemble)[i]->back().evaluate_ensemble_likelihood_ratios(index,
                                                                                                               inverse_incremental_temperature,
                                                                                                    conditioned_on_parameters);
  }
}
*/

void SequentialEnsembleKalmanWorker::subsample_weight(Ensemble* ensemble,
                                                      const Index* index,
                                                      double inverse_incremental_temperature)
{
  for (size_t i = 0; i < ensemble->size(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = (*ensemble)[i]->back().subsample_evaluate_ensemble_likelihood_ratios(index,
                                                                                                                         inverse_incremental_temperature);
  }
}

/*
void SequentialEnsembleKalmanWorker::subsample_weight(Ensemble* ensemble,
                                                      const Index* index,
                                                      double inverse_incremental_temperature,
                                                      const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < ensemble->size(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = (*ensemble)[i]->back().subsample_evaluate_ensemble_likelihood_ratios(index,
                                                                                                                                                 inverse_incremental_temperature,
                                                                                                                          conditioned_on_parameters);
  }
}
*/

arma::colvec SequentialEnsembleKalmanWorker::get_unnormalised_log_incremental_weights() const
{
  return this->log_unnormalised_incremental_weights;
}

/*
arma::colvec SequentialEnsembleKalmanWorker::get_unnormalised_log_incremental_weights() const
{
  return this->log_unnormalised_incremental_weights;
}
*/

void SequentialEnsembleKalmanWorker::specific_simulate(Ensemble* next_ensemble,
                                                       const Index* index)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    Particle* new_particle = next_ensemble->add_ensemble_member();
    this->the_enk->simulate_ensemble_member(local_rng,
                                            new_particle);
    //next_ensemble->push_back();
    
    //if (this->the_enk->likelihood_is_evaluated==true)
    //{
    //  next_ensemble->back()->back().evaluate_all_likelihoods(index);
    //}
  }
}

void SequentialEnsembleKalmanWorker::specific_simulate(Ensemble* next_ensemble,
                                                       const Index* index,
                                                       const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    Particle* new_particle = next_ensemble->add_ensemble_member();
    this->the_enk->simulate_ensemble_member(local_rng,
                                            new_particle,
                                            conditioned_on_parameters);
    
    //next_ensemble->push_back(this->the_enk->simulate_ensemble_member(local_rng,
    //                                                                 conditioned_on_parameters));
    
    //if (this->the_enk->likelihood_is_evaluated==true)
    //{
    //  next_ensemble->back()->back().evaluate_all_likelihoods(index,
    //                                                         conditioned_on_parameters);
    //}
  }
  
  //this->output = EnsembleKalman(this->particles);
}

/*
void SequentialEnsembleKalmanWorker::subsample_specific_simulate(EnsembleKalman* next_particles,
                                                      const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    next_particles->push_back(this->particle_simulator->subsample_simulate(local_rng,
                                                                           conditioned_on_parameters));
    
    if (this->the_enk->evaluate_log_proposal!=NULL)
    {
      (*next_particles)[i]->back().previous_target_evaluated = this->the_enk->evaluate_log_proposal((*next_particles)[i]->back().parameters);
    }
    
  }
  
  //this->output = EnsembleKalman(this->particles);
}

void SequentialEnsembleKalmanWorker::weight(EnsembleKalman &current_particles)
{
  
  //std::vector<LikelihoodEstimatorOutput*> inner_likelihood_estimator_outputs;
  //inner_likelihood_estimator_outputs.resize(this->get_number_of_particles());
  
  //likelihood_estimator_outputs.reserve(likelihood_estimators.size());
  //for (size_t j = 0; j < this->likelihood_estimators.size(); ++j)
  //{
  //  likelihood_estimator_outputs.push_back(inner_likelihood_estimator_outputs);
  //}
  
  
  
  // Set up likelihood estimators with all particles.
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator l = this->likelihood_estimator_outputs.begin();
  //     l!=this->likelihood_estimator_outputs.end();
  //     ++l)
  //{
  //  // Set up each estimator.
  //  // Not done yet.
  //}
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods() - current_particles[i]->back().previous_target_evaluated;
  }
  
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::iterator l_out = this->likelihood_estimator_outputs.begin();
  //     l_out!=likelihood_estimator_outputs.end();
  //     ++l_out)
  //{
  //  for (size_t l_in = 0;
  //       l_in!=this->get_number_of_particles();
  //        ++l_in)
  //   {
  //     (*l_out)[l_in]->estimate(this->particles[l_in].parameters);
  //     this->log_unnormalised_incremental_weights[l_in] = this->log_unnormalised_incremental_weights[l_in] + (*l_out)[l_in]->log_likelihood;
  //   }
  //}
  
}

void SequentialEnsembleKalmanWorker::pf_initial_weight(EnsembleKalman &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(&VectorSingleIndex(0)) - current_particles[i]->back().previous_target_evaluated;
  }
}

void SequentialEnsembleKalmanWorker::weight(EnsembleKalman &current_particles,
                                 const Parameters &conditioned_on_parameters)
{

  //std::vector<LikelihoodEstimatorOutput*> inner_likelihood_estimator_outputs;
  //inner_likelihood_estimator_outputs.resize(this->get_number_of_particles());
  
  //likelihood_estimator_outputs.reserve(likelihood_estimators.size());
  //for (size_t j = 0; j < this->likelihood_estimators.size(); ++j)
  //{
  //  likelihood_estimator_outputs.push_back(inner_likelihood_estimator_outputs);
  //}
  
  // Set up likelihood estimators with all particles.
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator l = this->likelihood_estimator_outputs.begin();
  //     l!=this->likelihood_estimator_outputs.end();
  //     ++l)
  //{
  //  // Set up each estimator.
  //  // Not done yet.
  //}
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
  }
  
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::iterator l_out = this->likelihood_estimator_outputs.begin();
  //     l_out!=likelihood_estimator_outputs.end();
  //     ++l_out)
  //{
  //  for (size_t l_in = 0;
  //       l_in!=this->get_number_of_particles();
  //        ++l_in)
  //   {
  //     (*l_out)[l_in]->estimate(this->particles[l_in].parameters);
  //     this->log_unnormalised_incremental_weights[l_in] = this->log_unnormalised_incremental_weights[l_in] + (*l_out)[l_in]->log_likelihood;
  //   }
  //}

}

void SequentialEnsembleKalmanWorker::pf_initial_weight(EnsembleKalman &current_particles,
                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
  }
}

void SequentialEnsembleKalmanWorker::subsample_weight(EnsembleKalman &current_particles,
                                           const Parameters &conditioned_on_parameters)
{
  
  //std::vector<LikelihoodEstimatorOutput*> inner_likelihood_estimator_outputs;
  //inner_likelihood_estimator_outputs.resize(this->get_number_of_particles());
  
  //likelihood_estimator_outputs.reserve(likelihood_estimators.size());
  //for (size_t j = 0; j < this->likelihood_estimators.size(); ++j)
  //{
  //  likelihood_estimator_outputs.push_back(inner_likelihood_estimator_outputs);
  //}
  
  
  
  // Set up likelihood estimators with all particles.
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator l = this->likelihood_estimator_outputs.begin();
  //     l!=this->likelihood_estimator_outputs.end();
  //     ++l)
  //{
  //  // Set up each estimator.
  //  // Not done yet.
  //}
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(conditioned_on_parameters) - current_particles[i]->back().subsample_previous_target_evaluated;
  }
  
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::iterator l_out = this->likelihood_estimator_outputs.begin();
  //     l_out!=likelihood_estimator_outputs.end();
  //     ++l_out)
  //{
  //  for (size_t l_in = 0;
  //       l_in!=this->get_number_of_particles();
  //        ++l_in)
  //   {
  //     (*l_out)[l_in]->estimate(this->particles[l_in].parameters);
  //     this->log_unnormalised_incremental_weights[l_in] = this->log_unnormalised_incremental_weights[l_in] + (*l_out)[l_in]->log_likelihood;
  //   }
  //}
  
}

void SequentialEnsembleKalmanWorker::smcfixed_weight(EnsembleKalman &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().evaluate_smcfixed_part_of_likelihoods();
  }
}

void SequentialEnsembleKalmanWorker::smcadaptive_given_smcfixed_weight(EnsembleKalman &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods() - current_particles[i]->back().previous_target_evaluated;
  }
}

void SequentialEnsembleKalmanWorker::smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods();
    
  }
}

void SequentialEnsembleKalmanWorker::marginal_weight(EnsembleKalman &current_particles,
                                          EnsembleKalman &previous_particles,
                                          ProposalKernel* proposal_kernel)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = previous_particles.normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                 previous_particles[j]->back());
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods() - log_sum_exp(terms);
    
  }
}

void SequentialEnsembleKalmanWorker::generic_weight(EnsembleKalman &current_particles,
                                         EnsembleKalman &previous_particles,
                                         ProposalKernel* proposal_kernel,
                                         ProposalKernel* L_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods()
    + L_kernel->evaluate_kernel(previous_particles[i]->back(),
                               current_particles[i]->back())
    - previous_particles[i]->back().target_evaluated
    - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                      previous_particles[i]->back());
    
  }
}

void SequentialEnsembleKalmanWorker::pf_weight(EnsembleKalman &current_particles,
                                    EnsembleKalman &previous_particles,
                                    ProposalKernel* proposal_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods();
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                       previous_particles[i]->back());
    
  }
}
*/

void SequentialEnsembleKalmanWorker::specific_move(Ensemble* next_particles,
                                                   Ensemble* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back()));
    
  }
}

/*
void SequentialEnsembleKalmanWorker::specific_move(Ensemble* next_particles,
                                                   Ensemble* current_particles,
                                                   const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back(),
                                                  conditioned_on_parameters));
    
  }
}
*/

void SequentialEnsembleKalmanWorker::subsample_specific_move(Ensemble* next_particles,
                                                             Ensemble* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back()));
    
  }
}

/*
void SequentialEnsembleKalmanWorker::subsample_specific_move(Ensemble* next_particles,
                                                             Ensemble* current_particles,
                                                             const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_ensemble_members());
  
  for (size_t i = 0; i < this->get_number_of_ensemble_members(); ++i)
  {
    next_particles->push_back(this->the_enk->move(local_rng,
                                                  (*current_particles)[i]->back(),
                                                  conditioned_on_parameters));
    
  }
}
*/

/*
void SequentialEnsembleKalmanWorker::smcfixed_weight(EnsembleKalman &current_particles,
                                          const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().evaluate_smcfixed_part_of_likelihoods(conditioned_on_parameters);
  }
}

void SequentialEnsembleKalmanWorker::smcadaptive_given_smcfixed_weight(EnsembleKalman &current_particles,
                                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
    
  }
}

void SequentialEnsembleKalmanWorker::smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_particles,
                                                                     const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(conditioned_on_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::marginal_weight(EnsembleKalman &current_particles,
                                          EnsembleKalman &previous_particles,
                                          ProposalKernel* proposal_kernel,
                                          const Parameters &conditioned_on_parameters)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = previous_particles.normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                 previous_particles[j]->back(),
                                                 conditioned_on_parameters);
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(conditioned_on_parameters) - log_sum_exp(terms);
    
  }
}

void SequentialEnsembleKalmanWorker::generic_weight(EnsembleKalman &current_particles,
                                         EnsembleKalman &previous_particles,
                                         ProposalKernel* proposal_kernel,
                                         ProposalKernel* L_kernel,
                                         const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(conditioned_on_parameters)
    + L_kernel->evaluate_kernel(previous_particles[i]->back(),
                               current_particles[i]->back(),
                               conditioned_on_parameters)
    - previous_particles[i]->back().target_evaluated
    - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                      previous_particles[i]->back(),
                                      conditioned_on_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::pf_weight(EnsembleKalman &current_particles,
                                    EnsembleKalman &previous_particles,
                                    ProposalKernel* proposal_kernel,
                                    const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(conditioned_on_parameters);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                                                       previous_particles[i]->back(),
                                                                                                                                       conditioned_on_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_specific_move(EnsembleKalman* next_particles,
                                                  const EnsembleKalman* current_particles,
                                                  const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    next_particles->push_back(this->the_enk->subsample_move(local_rng,
                                                            (*current_particles)[i]->back(),
                                                            conditioned_on_parameters));
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_smcfixed_weight(EnsembleKalman &current_particles,
                                          const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().subsample_evaluate_smcfixed_part_of_likelihoods(conditioned_on_parameters);
  }
}

void SequentialEnsembleKalmanWorker::subsample_smcadaptive_given_smcfixed_weight(EnsembleKalman &current_particles,
                                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(conditioned_on_parameters) - current_particles[i]->back().subsample_previous_target_evaluated;
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_smcadaptive_given_smcfixed_evaluate_target(EnsembleKalman &current_particles,
                                                                     const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(conditioned_on_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_marginal_weight(EnsembleKalman &current_particles,
                                          EnsembleKalman &previous_particles,
                                          ProposalKernel* proposal_kernel,
                                          const Parameters &conditioned_on_parameters)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = previous_particles.normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                 previous_particles[j]->back(),
                                                                                                 conditioned_on_parameters);
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(conditioned_on_parameters) - log_sum_exp(terms);
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_generic_weight(EnsembleKalman &current_particles,
                                         EnsembleKalman &previous_particles,
                                         ProposalKernel* proposal_kernel,
                                         ProposalKernel* L_kernel,
                                         const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(conditioned_on_parameters)
    + L_kernel->subsample_evaluate_kernel(previous_particles[i]->back(),
                                          current_particles[i]->back(),
                                          conditioned_on_parameters)
    - previous_particles[i]->back().subsample_target_evaluated
    - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                 previous_particles[i]->back(),
                                                 conditioned_on_parameters);
    
  }
}

void SequentialEnsembleKalmanWorker::subsample_pf_weight(EnsembleKalman &current_particles,
                                              EnsembleKalman &previous_particles,
                                              ProposalKernel* proposal_kernel,
                                              const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(conditioned_on_parameters);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                                                                                                       previous_particles[i]->back(),
                                                                                                                                       conditioned_on_parameters);
    
  }
}
*/

//void SequentialEnsembleKalmanWorker::specific_simulate_and_weight()
//{
//  RandomNumberGenerator local_rng(*this->get_rng());
//  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
//  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
//  {
//    this->particles[i] = (*this->particle_simulator)(local_rng);
//  }
  
  //std::vector<LikelihoodEstimatorOutput*> inner_likelihood_estimator_outputs;
  //inner_likelihood_estimator_outputs.resize(this->get_number_of_particles());
  
  //likelihood_estimator_outputs.reserve(likelihood_estimators.size());
  //for (size_t j = 0; j < this->likelihood_estimators.size(); ++j)
  //{
  //  likelihood_estimator_outputs.push_back(inner_likelihood_estimator_outputs);
  //}
  
//  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
//  {
//    this->log_unnormalised_incremental_weights[i] = 0.0;
//    for (std::vector<EvaluateLogDistributionPtr>::const_iterator p = this->prior_evaluates.begin();
//         p!=prior_evaluates.end();
//         ++p)
//    {
//      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + (*p)(this->particles[i].parameters);
//    }
    
//    for (std::vector<EvaluateLogDistributionPtr>::const_iterator p = this->proposal_evaluates.begin();
//         p!=proposal_evaluates.end();
//         ++p)
//    {
//      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - (*p)(this->particles[i].parameters);
//    }
    
    //size_t j = 0;
    //std::vector<LikelihoodEstimatorOutput*> outputs;
    //outputs.reserve(likelihood_estimators.size());
    //for (std::vector<LikelihoodEstimator*>::const_iterator l = this->likelihood_estimators.begin();
    //     l!=likelihood_estimators.end();
    //     ++l)
    //{
    //  outputs.push_back((*l)->initial_simulate(this->particles[i].parameters));
    //}
    
//  }
  
  // Set up likelihood estimators with all particles.
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::const_iterator l = this->likelihood_estimator_outputs.begin();
  //     l!=this->likelihood_estimator_outputs.end();
  //     ++l)
  //{
  //  // Set up each estimator.
  //  // Not done yet.
  //}
  
//  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
//  {
//    this->particles[i].estimate_likelihoods();
//    this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + this->particles[i].get_log_likelihood();
//  }
  
  //for (std::vector< std::vector<LikelihoodEstimatorOutput*> >::iterator l_out = this->likelihood_estimator_outputs.begin();
  //     l_out!=likelihood_estimator_outputs.end();
  //     ++l_out)
  //{
  //  for (size_t l_in = 0;
  //       l_in!=this->get_number_of_particles();
  //        ++l_in)
  //   {
  //     (*l_out)[l_in]->estimate(this->particles[l_in].parameters);
  //     this->log_unnormalised_incremental_weights[l_in] = this->log_unnormalised_incremental_weights[l_in] + (*l_out)[l_in]->log_likelihood;
  //   }
  //}

//}
