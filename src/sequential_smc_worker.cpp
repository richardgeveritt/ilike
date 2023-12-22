#include "sequential_smc_worker.h"
#include "particle_simulator.h"
#include "smc.h"
#include "likelihood_estimator_output.h"
#include "utils.h"
#include "move_output.h"
#include "index.h"
#include "vector_single_index.h"
#include "factors.h"
#include "factor_variables.h"

//Default constructor.
SequentialSMCWorker::SequentialSMCWorker()
  :SMCWorker()
{
}

SequentialSMCWorker::SequentialSMCWorker(SMC* the_smc_in)
  :SMCWorker(the_smc_in)
{
  //this->particles = NULL;
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
  
  //this->particles.clear();

  SMCWorker::operator=(another);
  this->make_copy(another);
}

SMCWorker* SequentialSMCWorker::duplicate() const
{
  return( new SequentialSMCWorker(*this));
}

/*
Particles SequentialSMCWorker::simulated_particles() const
{
  return this->particles;
}

Particles& SequentialSMCWorker::simulated_particles()
{
  return this->particles;
}
*/

void SequentialSMCWorker::make_copy(const SequentialSMCWorker &another)
{
  //this->particles = another.particles;
  this->log_unnormalised_incremental_weights = another.log_unnormalised_incremental_weights;
}

arma::colvec SequentialSMCWorker::get_unnormalised_log_incremental_weights() const
{
  return this->log_unnormalised_incremental_weights;
}

void SequentialSMCWorker::specific_simulate(Particles* next_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    //Particle* new_particle = next_particles->add_particle();
    //this->the_smc->particle_simulator->simulate(local_rng,
    //                                            new_particle,
    //                                            this->the_smc->factors);
    
    next_particles->push_back(this->the_smc->particle_simulator->simulate(local_rng,
                                                                          this->the_smc->factors,
                                                                          &this->the_smc->proposals_to_transform_for,
                                                                          &this->the_smc->proposals_to_find_gradient_for,
                                                                          this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
    }
    
  }
  
  /*
  if (this->the_smc->transform==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      Particle* new_particle = next_particles->add_particle();
      this->the_smc->particle_simulator->simulate(local_rng,
                                                  new_particle,
                                                  this->the_smc->factors);
      
      //next_particles->push_back(new_p);
      
      if (this->the_smc->proposal_is_evaluated==true)
      {
        next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
      }
      
    }
  }
  else
  {
    if (this->the_smc->store_transformed==false)
    {
      Rcpp::stop("SequentialSMCWorker::specific_simulate - need to store transformed auxiliary variables if transform is specified.");
    }
    else
    {
      for (size_t i = 0; i < this->get_number_of_particles(); ++i)
      {
        Particle* new_particle = next_particles->add_particle();
        this->the_smc->particle_simulator->simulate_and_transform(local_rng,
                                                                  new_particle,
                                                                  this->the_smc->factors,
                                                                  this->the_smc->transform,
                                                                  this->the_smc->store_raw);
        
        if (this->the_smc->proposal_is_evaluated==true) // if evaluated, proposal must be on transformed space
        {
          next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
        }

      }
    }
  }
  */

  //this->output = Particles(this->particles);
  
  /*
   for (size_t i = 0; i < this->get_number_of_particles(); ++i)
   {
   Parameters merged_parameters = this->particles[i].parameters.merge(conditioned_on_parameters);
   this->log_unnormalised_incremental_weights[i] = 0.0;
   if (this->evaluate_log_prior!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + this->evaluate_log_prior(merged_parameters);
   if (this->evaluate_log_proposal!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - this->evaluate_log_proposal(merged_parameters);
   
   //size_t j = 0;
   //std::vector<LikelihoodEstimatorOutput*> outputs;
   //outputs.reserve(likelihood_estimators.size());
   //for (std::vector<LikelihoodEstimator*>::const_iterator l = this->likelihood_estimators.begin();
   //     l!=likelihood_estimators.end();
   //     ++l)
   //{
   //  outputs.push_back((*l)->initial_simulate(this->particles[i].parameters));
   //}
   }
   */
}

void SequentialSMCWorker::specific_simulate(Particles* next_particles,
                                            const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    //Particle* new_particle = next_particles->add_particle();
    //this->the_smc->particle_simulator->simulate(local_rng,
    //                                            new_particle,
    //                                            this->the_smc->factors,
    //                                            conditioned_on_parameters);
    
    next_particles->push_back(this->the_smc->particle_simulator->simulate(local_rng,
                                                                          this->the_smc->factors,
                                                                          &this->the_smc->proposals_to_transform_for,
                                                                          &this->the_smc->proposals_to_find_gradient_for,
                                                                          conditioned_on_parameters,
                                                                          this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back());
    }
    
  }
  
  /*
  if (this->the_smc->transform==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      Particle* new_particle = next_particles->add_particle();
      this->the_smc->particle_simulator->simulate(local_rng,
                                                  new_particle,
                                                  this->the_smc->factors,
                                                  conditioned_on_parameters);
      
      //next_particles->push_back(this->the_smc->particle_simulator->simulate(local_rng,
      //                                                                      this->the_smc->factors,
      //                                                                      conditioned_on_parameters));
      
      if (this->the_smc->proposal_is_evaluated==true)
      {
        next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back(),conditioned_on_parameters);
      }
      
    }
  }
  else
  {
    if (this->the_smc->store_transformed==false)
    {
      Rcpp::stop("SequentialSMCWorker::specific_simulate - need to store transformed auxiliary variables if transform is specified.");
    }
    else
    {
      for (size_t i = 0; i < this->get_number_of_particles(); ++i)
      {
        Particle* new_particle = next_particles->add_particle();
        this->the_smc->particle_simulator->simulate_and_transform(local_rng,
                                                                  new_particle,
                                                                  this->the_smc->factors,
                                                                  this->the_smc->transform,
                                                                  this->the_smc->store_raw,
                                                                  conditioned_on_parameters);
        
        if (this->the_smc->proposal_is_evaluated==true) // if evaluated, proposal must be on transformed space
        {
          next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->evaluate(next_particles->back()->back(),conditioned_on_parameters);
        }
        
      }
    }
  }
  */
  
  //this->output = Particles(this->particles);
  
  /*
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    Parameters merged_parameters = this->particles[i].parameters.merge(conditioned_on_parameters);
    this->log_unnormalised_incremental_weights[i] = 0.0;
    if (this->evaluate_log_prior!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + this->evaluate_log_prior(merged_parameters);
    if (this->evaluate_log_proposal!=NULL)
    this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - this->evaluate_log_proposal(merged_parameters);
    
    //size_t j = 0;
    //std::vector<LikelihoodEstimatorOutput*> outputs;
    //outputs.reserve(likelihood_estimators.size());
    //for (std::vector<LikelihoodEstimator*>::const_iterator l = this->likelihood_estimators.begin();
    //     l!=likelihood_estimators.end();
    //     ++l)
    //{
    //  outputs.push_back((*l)->initial_simulate(this->particles[i].parameters));
    //}
  }
  */
}

void SequentialSMCWorker::subsample_specific_simulate(Particles* next_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    //Particle* new_particle = next_particles->add_particle();
    //this->the_smc->particle_simulator->subsample_simulate(local_rng,
    //                                                      new_particle,
    //                                                      this->the_smc->factors,
    //                                                      conditioned_on_parameters);
    
    next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                    this->the_smc->factors,
                                                                                    &this->the_smc->proposals_to_transform_for,
                                                                                    &this->the_smc->proposals_to_find_gradient_for,
                                                                                    this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back());
    }
    
  }
  
  /*
   if (this->the_smc->transform==NULL)
   {
   for (size_t i = 0; i < this->get_number_of_particles(); ++i)
   {
   Particle* new_particle = next_particles->add_particle();
   this->the_smc->particle_simulator->subsample_simulate(local_rng,
   new_particle,
   this->the_smc->factors,
   conditioned_on_parameters);
   
   //next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
   //                                                                                this->the_smc->factors,
   //                                                                                conditioned_on_parameters));
   
   if (this->the_smc->proposal_is_evaluated==true)
   {
   next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back(),conditioned_on_parameters);
   }
   
   }
   }
   else
   {
   if (this->the_smc->store_transformed==false)
   {
   Rcpp::stop("SequentialSMCWorker::specific_simulate - need to store transformed auxiliary variables if transform is specified.");
   }
   else
   {
   for (size_t i = 0; i < this->get_number_of_particles(); ++i)
   {
   Particle* new_particle = next_particles->add_particle();
   this->the_smc->particle_simulator->subsample_simulate_and_transform(local_rng,
   new_particle,
   this->the_smc->factors,
   this->the_smc->transform,
   this->the_smc->store_raw,
   conditioned_on_parameters);
   
   if (this->the_smc->proposal_is_evaluated==true) // if evaluated, proposal must be on transformed space
   {
   next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back(),conditioned_on_parameters);
   }
   
   }
   }
   }
   */
  
  //this->output = Particles(this->particles);
  
  /*
   for (size_t i = 0; i < this->get_number_of_particles(); ++i)
   {
   Parameters merged_parameters = this->particles[i].parameters.merge(conditioned_on_parameters);
   this->log_unnormalised_incremental_weights[i] = 0.0;
   if (this->evaluate_log_prior!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + this->evaluate_log_prior(merged_parameters);
   if (this->evaluate_log_proposal!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - this->evaluate_log_proposal(merged_parameters);
   
   //size_t j = 0;
   //std::vector<LikelihoodEstimatorOutput*> outputs;
   //outputs.reserve(likelihood_estimators.size());
   //for (std::vector<LikelihoodEstimator*>::const_iterator l = this->likelihood_estimators.begin();
   //     l!=likelihood_estimators.end();
   //     ++l)
   //{
   //  outputs.push_back((*l)->initial_simulate(this->particles[i].parameters));
   //}
   }
   */
}

void SequentialSMCWorker::subsample_specific_simulate(Particles* next_particles,
                                                      const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    //Particle* new_particle = next_particles->add_particle();
    //this->the_smc->particle_simulator->subsample_simulate(local_rng,
    //                                                      new_particle,
    //                                                      this->the_smc->factors,
    //                                                      conditioned_on_parameters);
    
    next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                                                    this->the_smc->factors,
                                                                                    &this->the_smc->proposals_to_transform_for,
                                                                                    &this->the_smc->proposals_to_find_gradient_for,
                                                                                    conditioned_on_parameters,
                                                                                    this->the_smc->sequencer.schedule_parameters));
    
    if (this->the_smc->proposal_is_evaluated==true)
    {
      next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back());
    }
    
  }
  
  /*
  if (this->the_smc->transform==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      Particle* new_particle = next_particles->add_particle();
      this->the_smc->particle_simulator->subsample_simulate(local_rng,
                                                            new_particle,
                                                            this->the_smc->factors,
                                                            conditioned_on_parameters);
      
      //next_particles->push_back(this->the_smc->particle_simulator->subsample_simulate(local_rng,
      //                                                                                this->the_smc->factors,
      //                                                                                conditioned_on_parameters));
      
      if (this->the_smc->proposal_is_evaluated==true)
      {
        next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back(),conditioned_on_parameters);
      }
      
    }
  }
  else
  {
    if (this->the_smc->store_transformed==false)
    {
      Rcpp::stop("SequentialSMCWorker::specific_simulate - need to store transformed auxiliary variables if transform is specified.");
    }
    else
    {
      for (size_t i = 0; i < this->get_number_of_particles(); ++i)
      {
        Particle* new_particle = next_particles->add_particle();
        this->the_smc->particle_simulator->subsample_simulate_and_transform(local_rng,
                                                                            new_particle,
                                                                            this->the_smc->factors,
                                                                            this->the_smc->transform,
                                                                            this->the_smc->store_raw,
                                                                            conditioned_on_parameters);
        
        if (this->the_smc->proposal_is_evaluated==true) // if evaluated, proposal must be on transformed space
        {
          next_particles->back()->back().previous_target_evaluated = this->the_smc->particle_simulator->subsample_evaluate(next_particles->back()->back(),conditioned_on_parameters);
        }
        
      }
    }
  }
  */
  
  //this->output = Particles(this->particles);
  
  /*
   for (size_t i = 0; i < this->get_number_of_particles(); ++i)
   {
   Parameters merged_parameters = this->particles[i].parameters.merge(conditioned_on_parameters);
   this->log_unnormalised_incremental_weights[i] = 0.0;
   if (this->evaluate_log_prior!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] + this->evaluate_log_prior(merged_parameters);
   if (this->evaluate_log_proposal!=NULL)
   this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - this->evaluate_log_proposal(merged_parameters);
   
   //size_t j = 0;
   //std::vector<LikelihoodEstimatorOutput*> outputs;
   //outputs.reserve(likelihood_estimators.size());
   //for (std::vector<LikelihoodEstimator*>::const_iterator l = this->likelihood_estimators.begin();
   //     l!=likelihood_estimators.end();
   //     ++l)
   //{
   //  outputs.push_back((*l)->initial_simulate(this->particles[i].parameters));
   //}
   }
   */
}

void SequentialSMCWorker::weight(const Index* index,
                                 Particles &current_particles)
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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index) - current_particles[i]->back().previous_target_evaluated;
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

void SequentialSMCWorker::pf_initial_weight(Particles &current_particles)
{
  VectorSingleIndex index(0);
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(&index) - current_particles[i]->back().previous_target_evaluated;
  }
}

/*
void SequentialSMCWorker::weight(const Index* index,
                                 Particles &current_particles,
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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index,conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
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
*/

/*
void SequentialSMCWorker::pf_initial_weight(Particles &current_particles,
                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    VectorSingleIndex index(0);
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(&index,
                                                                                                      conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
  }
}
*/

void SequentialSMCWorker::subsample_weight(const Index* index,
                                           Particles &current_particles)
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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index) - current_particles[i]->back().subsample_previous_target_evaluated;
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

/*
void SequentialSMCWorker::subsample_weight(const Index* index,
                                           Particles &current_particles,
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
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index,conditioned_on_parameters) - current_particles[i]->back().subsample_previous_target_evaluated;
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
*/

void SequentialSMCWorker::subsample_pf_initial_weight(Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    VectorSingleIndex index(0);
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(&index) - current_particles[i]->back().previous_target_evaluated;
  }
}

/*
void SequentialSMCWorker::subsample_pf_initial_weight(Particles &current_particles,
                                                      const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    VectorSingleIndex index(0);
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(&index,
                                                                                                                conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
  }
}
*/

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
    //double prev = current_particles[i]->back().previous_target_evaluated;
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
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    next_particles->push_back(this->the_smc->move(local_rng,
                                                  (*current_particles)[current_particles->ancestor_variables[i]]->back()));
    //Particle* new_particle = next_particles->add_particle();
    //this->the_smc->move(local_rng,
    //                    new_particle,
    //                    (*current_particles)[current_particles->ancestor_variables[i]]->back());
    
  }
  
  /*
  if (this->the_smc->sequencer_parameters==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->move(local_rng,
                                                    (*current_particles)[current_particles->ancestor_variables[i]]->back()));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->move(local_rng,
      //                    new_particle,
      //                    (*current_particles)[current_particles->ancestor_variables[i]]->back());
      
    }
  }
  else
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->move(local_rng,
                                                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                    *this->the_smc->sequencer_parameters));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->move(local_rng,
      //                    next_particle,
      //                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                    *this->the_smc->sequencer_parameters);
      
    }
  }
  */
}

/*
void SequentialSMCWorker::smcfixed_weight(const Index* index,
                                          Particles &current_particles,
                                          const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().evaluate_smcfixed_part_of_likelihoods(index,
                                                                       conditioned_on_parameters);
  }
}
*/

/*
void SequentialSMCWorker::smcadaptive_given_smcfixed_weight(const Index* index,
                                                            Particles &current_particles,
                                                            const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,conditioned_on_parameters) - current_particles[i]->back().previous_target_evaluated;
    
  }
}
*/

/*
void SequentialSMCWorker::smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                     Particles &current_particles,
                                                                     const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                                                                      conditioned_on_parameters);
    
  }
}
*/

/*
void SequentialSMCWorker::marginal_weight(const Index* index,
                                          Particles &current_particles,
                                          Particles &previous_particles,
                                          ProposalKernel* proposal_kernel,
                                          const Parameters &conditioned_on_parameters)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = current_particles.previous_normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                 previous_particles[previous_particles.ancestor_variables[j]]->back(),
                                                 conditioned_on_parameters);
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index,
                                                                                                      conditioned_on_parameters) - log_sum_exp(terms);
    
  }
}
*/

/*
void SequentialSMCWorker::generic_weight(const Index* index,
                                         Particles &current_particles,
                                         Particles &previous_particles,
                                         ProposalKernel* proposal_kernel,
                                         ProposalKernel* L_kernel,
                                         const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index,
                                                                                                      conditioned_on_parameters)
    + L_kernel->evaluate_kernel(previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                current_particles[i]->back(),
                                conditioned_on_parameters)
    - previous_particles[previous_particles.ancestor_variables[i]]->back().target_evaluated
    - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                       previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                       conditioned_on_parameters);
    
  }
}
*/

/*
void SequentialSMCWorker::pf_weight(const Index* index,
                                    Particles &current_particles,
                                    Particles &previous_particles,
                                    ProposalKernel* proposal_kernel,
                                    const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().evaluate_likelihoods(index,
                                                                                                      conditioned_on_parameters);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                                                       previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                                                                                                                       conditioned_on_parameters);
    
  }
}
*/

/*
void SequentialSMCWorker::specific_move(Particles* next_particles,
                                        const Particles* current_particles,
                                        const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  if (this->the_smc->sequencer_parameters==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->move(local_rng,
                                                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                    conditioned_on_parameters));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->move(local_rng,
      //                    new_particle,
      //                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                    conditioned_on_parameters);
      
    }
  }
  else
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->move(local_rng,
                                                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                    conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters)));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->move(local_rng,
      //                    new_particle,
      //                    (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                    conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters));
      
    }
  }
}
*/

void SequentialSMCWorker::subsample_specific_move(Particles* next_particles,
                                                  const Particles* current_particles)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                            (*current_particles)[current_particles->ancestor_variables[i]]->back()));
  }
  
  /*
  if (this->the_smc->sequencer_parameters==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                              (*current_particles)[current_particles->ancestor_variables[i]]->back()));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->subsample_move(local_rng,
      //                              new_particle,
      //                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                              conditioned_on_parameters));
      
    }
  }
  else
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                              conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters)));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->subsample_move(local_rng,
      //                              new_particle,
      //                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                              conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters)));
      
    }
  }
  */
}

/*
void SequentialSMCWorker::subsample_specific_move(Particles* next_particles,
                                                  const Particles* current_particles,
                                                  const Parameters &conditioned_on_parameters)
{
  RandomNumberGenerator local_rng(*this->get_rng());
  local_rng.seed(this->get_seed(),this->get_number_of_particles());
  
  if (this->the_smc->sequencer_parameters==NULL)
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                              conditioned_on_parameters));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->subsample_move(local_rng,
      //                              new_particle,
      //                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                              conditioned_on_parameters));
      
    }
  }
  else
  {
    for (size_t i = 0; i < this->get_number_of_particles(); ++i)
    {
      next_particles->push_back(this->the_smc->subsample_move(local_rng,
                                                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
                                                              conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters)));
      //Particle* new_particle = next_particles->add_particle();
      //this->the_smc->subsample_move(local_rng,
      //                              new_particle,
      //                              (*current_particles)[current_particles->ancestor_variables[i]]->back(),
      //                              conditioned_on_parameters.merge(*this->the_smc->sequencer_parameters)));
      
    }
  }
}
*/

void SequentialSMCWorker::subsample_smcfixed_weight(const Index* index,
                                                    Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().subsample_evaluate_smcfixed_part_of_likelihoods(index);
  }
}

/*
void SequentialSMCWorker::subsample_smcfixed_weight(const Index* index,
                                                    Particles &current_particles,
                                                    const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    current_particles[i]->back().subsample_evaluate_smcfixed_part_of_likelihoods(index,
                                                                                 conditioned_on_parameters);
  }
}
*/

void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                                      Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index) - current_particles[i]->back().subsample_previous_target_evaluated;
    
  }
}

/*
void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                                      Particles &current_particles,
                                                                      const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                                                                                conditioned_on_parameters) - current_particles[i]->back().subsample_previous_target_evaluated;
    
  }
}
*/

void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                               Particles &current_particles)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);
    
  }
}

/*
void SequentialSMCWorker::subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                               Particles &current_particles,
                                                                               const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(index,
                                                                                                                                                conditioned_on_parameters);
    
  }
}
*/

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

/*
void SequentialSMCWorker::subsample_marginal_weight(const Index* index,
                                                    Particles &current_particles,
                                                    Particles &previous_particles,
                                                    ProposalKernel* proposal_kernel,
                                                    const Parameters &conditioned_on_parameters)
{
  arma::colvec terms(this->get_number_of_particles());
  
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    // If proposal is Gaussian, might be able to use the fast Gauss transform in low-dim
    for (size_t j = 0; j < this->get_number_of_particles(); ++j)
    {
      terms[j] = current_particles.previous_normalised_log_weights[j] + proposal_kernel->evaluate_kernel(current_particles[i]->back(),
                                                                                                 previous_particles[previous_particles.ancestor_variables[j]]->back(),
                                                                                                 conditioned_on_parameters);
    }
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index,
                                                                                                                conditioned_on_parameters) - log_sum_exp(terms);
    
  }
}
*/

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

/*
void SequentialSMCWorker::subsample_generic_weight(const Index* index,
                                                   Particles &current_particles,
                                                   Particles &previous_particles,
                                                   ProposalKernel* proposal_kernel,
                                                   ProposalKernel* L_kernel,
                                                   const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index,
                                                                                                                conditioned_on_parameters)
    + L_kernel->subsample_evaluate_kernel(previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                          current_particles[i]->back(),
                                          conditioned_on_parameters)
    - previous_particles[previous_particles.ancestor_variables[i]]->back().subsample_target_evaluated
    - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                 previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                                 conditioned_on_parameters);
    
  }
}
*/

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

/*
void SequentialSMCWorker::subsample_pf_weight(const Index* index,
                                              Particles &current_particles,
                                              Particles &previous_particles,
                                              ProposalKernel* proposal_kernel,
                                              const Parameters &conditioned_on_parameters)
{
  for (size_t i = 0; i < this->get_number_of_particles(); ++i)
  {
    
    this->log_unnormalised_incremental_weights[i] = current_particles[i]->back().subsample_evaluate_likelihoods(index,
                                                                                                                conditioned_on_parameters);
    if (proposal_kernel!=NULL)
      this->log_unnormalised_incremental_weights[i] = this->log_unnormalised_incremental_weights[i] - proposal_kernel->subsample_evaluate_kernel(current_particles[i]->back(),
                                                                                                                                       previous_particles[previous_particles.ancestor_variables[i]]->back(),
                                                                                                                                       conditioned_on_parameters);
    
  }
}
*/

//void SequentialSMCWorker::specific_simulate_and_weight()
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
