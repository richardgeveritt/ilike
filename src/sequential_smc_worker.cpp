#include "sequential_smc_worker.h"
#include "particle_simulator.h"
#include "smc.h"
#include "likelihood_estimator_output.h"
#include "utils.h"
#include "move_output.h"
#include "index.h"
#include "vector_index.h"
#include "factors.h"
#include "factor_variables.h"

namespace ilike
{
//Default constructor.
SequentialSMCWorker::SequentialSMCWorker()
:SMCWorker()
{
}

SequentialSMCWorker::SequentialSMCWorker(SMC* the_smc_in)
:SMCWorker(the_smc_in)
{
  this->log_unnormalised_incremental_weights = std::vector<double>(this->get_number_of_particles());
}

//Copy constructor for the SequentialSMCWorker class.
SequentialSMCWorker::SequentialSMCWorker(const SequentialSMCWorker &another)
:SMCWorker(another)
{
  this->make_copy(another);
}

//Destructor for the SequentialSMCWorker class.
SequentialSMCWorker::~SequentialSMCWorker()
{
}

void SequentialSMCWorker::operator=(const SequentialSMCWorker &another)
{
  if(this == &another){ //if a==a
    return;
  }

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* SequentialSMCWorker::duplicate() const
{
  return( new SequentialSMCWorker(*this));
}

void SequentialSMCWorker::make_copy(const SequentialSMCWorker &another)
{
  this->log_unnormalised_incremental_weights = another.log_unnormalised_incremental_weights;
}

arma::colvec SequentialSMCWorker::get_unnormalised_log_incremental_weights() const
{
  return this->log_unnormalised_incremental_weights;
}

void SequentialSMCWorker::set_unnormalised_log_incremental_weights(const arma::colvec &weights)
{
  this->log_unnormalised_incremental_weights = weights;
}

void SequentialSMCWorker::specific_simulate(Particles* next_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    next_particles->push_back(this->the_smc->particle_simulator->simulate(local_rng,
                                                                          this->the_smc->factors.get(),
                                                                          &this->the_smc->proposals_to_transform_for,
                                                                          &this->the_smc->proposals_to_find_gradient_for,
                                                                          this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
    }
    
  }

}

void SequentialSMCWorker::specific_simulate(Particles* next_particles,
                                            const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    next_particles->push_back(this->the_smc->particle_simulator->simulate(local_rng,
                                                                          this->the_smc->factors.get(),
                                                                          &this->the_smc->proposals_to_transform_for,
                                                                          &this->the_smc->proposals_to_find_gradient_for,
                                                                          conditioned_on_parameters,
                                                                          this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
    }
    
  }

}

void SequentialSMCWorker::subsample_specific_simulate(Particles* next_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                    this->the_smc->factors.get(),
                                                                                    &this->the_smc->proposals_to_transform_for,
                                                                                    &this->the_smc->proposals_to_find_gradient_for,
                                                                                    this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back());
    }
    
  }

}

void SequentialSMCWorker::subsample_specific_simulate(Particles* next_particles,
                                                      const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                    this->the_smc->factors.get(),
                                                                                    &this->the_smc->proposals_to_transform_for,
                                                                                    &this->the_smc->proposals_to_find_gradient_for,
                                                                                    conditioned_on_parameters,
                                                                                    this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back());
    }
    
  }

}

void SequentialSMCWorker::weight(const Index* index,
                                 Particles &current_particles)
{

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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index) - current_particles[i]->back().previous_target_evaluated;
  }

}

void SequentialSMCWorker::pf_initial_weight(Particles &current_particles)
{
  VectorIndex index(0);
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(&index) - current_particles[i]->back().previous_target_evaluated;
  }
}

void SequentialSMCWorker::subsample_weight(const Index* index,
                                           Particles &current_particles)
{

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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index) - current_particles[i]->back().subsample_previous_target_evaluated;
  }

}

void SequentialSMCWorker::subsample_pf_initial_weight(Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    VectorIndex index(0);
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(&index) - current_particles[i]->back().previous_target_evaluated;
  }
}

void SequentialSMCWorker::smcfixed_weight(const Index* index,
                                          Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().evaluate_smcfixed_part_of_likelihoods(index);
  }
}

void SequentialSMCWorker::smcadaptive_given_smcfixed_weight(const Index* index,
                                                            Particles &current_particles)
{
  size_t n = this->get_number_of_particles();
  for (size_t i = 0; i < n; ++i)
  {
    double a = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    double b = current_particles[i]->back().previous_target_evaluated;
    
    if (a==-arma::datum::inf)
    {
      this->log_unnormalised_incremental_weights[i] = -arma::datum::inf;
    }
    else
    {
      this->log_unnormalised_incremental_weights[i] = a - b;
    }
  }
}

void SequentialSMCWorker::smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                     Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    
  }
}

void SequentialSMCWorker::marginal_weight(const Index* index,
                                          Particles &current_particles,
                                          Particles &previous_particles,
                                          ProposalKernel* proposal_kernel)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = current_particles.previous_normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                         previous_particles[previous_particles.ancestor_variables[j]]->back());
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index) - log_sum_exp(terms);
    
  }
}

void SequentialSMCWorker::generic_weight(const Index* index,
                                         Particles &current_particles,
                                         Particles &previous_particles,
                                         ProposalKernel* proposal_kernel,
                                         ProposalKernel* L_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index)
    + L_kernel->evaluate_kernel(previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                current_particles[i]->back())
    - previous_particles[previous_particles.ancestor_variables[i]]->back().target_evaluated
    - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                       previous_particles[previous_particles.ancestor_variables[i]]->back());
    
  }
}

void SequentialSMCWorker::pf_weight(const Index* index,
                                    Particles &current_particles,
                                    Particles &previous_particles,
                                    ProposalKernel* proposal_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                                                       previous_particles[previous_particles.ancestor_variables[i]]->back());
    
  }
}

void SequentialSMCWorker::specific_move(Particles* next_particles,
                                        const Particles* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    MoveOutput* spare = next_particles->take_spare();
    if (spare!=NULL)
      next_particles->push_back(this->the_smc->move_reusing(local_rng,
                                                            (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                            spare));
    else
      next_particles->push_back(this->the_smc->move(local_rng,
                                                    (*current_particles)[current_particles->ancestor_variables[i]]->back()));
    next_particles->back()->back().previous_self = &(*current_particles)[current_particles->ancestor_variables[i]]->back();
    
  }
  
}

void SequentialSMCWorker::subsample_specific_move(Particles* next_particles,
                                                  const Particles* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                            (*current_particles)[current_particles->ancestor_variables[i]]->back()));
    next_particles->back()->back().previous_self = &(*current_particles)[current_particles->ancestor_variables[i]]->back();
  }
  
}

void SequentialSMCWorker::subsample_smcfixed_weight(const Index* index,
                                                    Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().subsample_evaluate_smcfixed_part_of_likelihoods(index);
  }
}

void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                                      Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index) - current_particles[i]->back().subsample_previous_target_evaluated;
    
  }
}

void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                               Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    
  }
}

void SequentialSMCWorker::subsample_marginal_weight(const Index* index,
                                                    Particles &current_particles,
                                                    Particles &previous_particles,
                                                    ProposalKernel* proposal_kernel)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = current_particles.previous_normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                         previous_particles[previous_particles.ancestor_variables[j]]->back());
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index) - log_sum_exp(terms);
    
  }
}

void SequentialSMCWorker::subsample_generic_weight(const Index* index,
                                                   Particles &current_particles,
                                                   Particles &previous_particles,
                                                   ProposalKernel* proposal_kernel,
                                                   ProposalKernel* L_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index)
    + L_kernel->subsample_evaluate_kernel(previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                          current_particles[i]->back())
    - previous_particles[previous_particles.ancestor_variables[i]]->back().subsample_target_evaluated
    - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                 previous_particles[previous_particles.ancestor_variables[i]]->back());
    
  }
}

void SequentialSMCWorker::subsample_pf_weight(const Index* index,
                                              Particles &current_particles,
                                              Particles &previous_particles,
                                              ProposalKernel* proposal_kernel)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                                                                                                                 previous_particles[previous_particles.ancestor_variables[i]]->back());
    
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
