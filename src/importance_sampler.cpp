#include "importance_sampler.h"
#include "smc_worker.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "exact_likelihood_estimator.h"
#include "parameter_particle_simulator.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "independent_proposal_kernel.h"
#include "custom_independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_single_index.h"

ImportanceSampler::ImportanceSampler()
  :SMC()
{
  this->index = NULL;
}

ImportanceSampler::ImportanceSampler(RandomNumberGenerator* rng_in,
                                     size_t* seed_in,
                                     Data* data_in,
                                     size_t number_of_particles_in,
                                     EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                     SimulateIndependentProposalPtr simulate_prior_in,
                                     bool smcfixed_flag_in,
                                     bool parallel_in,
                                     size_t grain_size_in)
  :SMC(rng_in, seed_in, data_in, number_of_particles_in, 1, 0, double(number_of_particles_in), false, smcfixed_flag_in, true)
{
  this->index = NULL;
  
  std::vector<size_t> indices;
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  likelihood_estimators.reserve(1);
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               evaluate_log_likelihood_in,
                                                               true));
  indices.push_back(0);
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_prior_in);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);

  if (parallel_in==TRUE)
  {
      //this->the_worker = new RcppParallelSMCWorker(this,
                                                //this->model_and_algorithm.particle_simulator,
                                                //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }

}

ImportanceSampler::ImportanceSampler(RandomNumberGenerator* rng_in,
                                     size_t* seed_in,
                                     Data* data_in,
                                     size_t number_of_particles_in,
                                     EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                     EvaluateLogDistributionPtr evaluate_log_prior_in,
                                     SimulateIndependentProposalPtr simulate_proposal_in,
                                     EvaluateLogDistributionPtr evaluate_log_proposal_in,
                                     bool smcfixed_flag_in,
                                     bool parallel_in,
                                     size_t grain_size_in)
  :SMC(rng_in, seed_in, data_in, number_of_particles_in, 1, 0, double(number_of_particles_in), true, smcfixed_flag_in, true)
{
  this->index = NULL;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  likelihood_estimators.reserve(1);
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               evaluate_log_prior_in,
                                                               evaluate_log_likelihood_in,
                                                               true));
  indices.push_back(0);
  this->index = new VectorSingleIndex(indices);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators,
  //                                                                              "u");

  if (parallel_in==TRUE)
  {
      //this->the_worker = new RcppParallelSMCWorker(this,
                                                //this->model_and_algorithm.particle_simulator,
                                                //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }

}

ImportanceSampler::ImportanceSampler(RandomNumberGenerator* rng_in,
                                     size_t* seed_in,
                                     Data* data_in,
                                     size_t number_of_particles_in,
                                     EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                     EvaluateLogDistributionPtr evaluate_log_prior_in,
                                     IndependentProposalKernel* proposal_in,
                                     bool smcfixed_flag_in,
                                     bool parallel_in,
                                     size_t grain_size_in)
:SMC(rng_in, seed_in, data_in, number_of_particles_in, 1, 0, double(number_of_particles_in), true, smcfixed_flag_in, true)
{
  this->index = NULL;
  
  std::vector<size_t> indices;
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  likelihood_estimators.reserve(1);
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               evaluate_log_prior_in,
                                                               evaluate_log_likelihood_in,
                                                               true));
  
  indices.push_back(0);
  this->index = new VectorSingleIndex(indices);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators);
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators,
  //                                                                              "u");
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,
    //this->model_and_algorithm.particle_simulator,
    //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
}

//Copy constructor for the ImportanceSampler class.
ImportanceSampler::ImportanceSampler(const ImportanceSampler &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the ImportanceSampler class.
ImportanceSampler::~ImportanceSampler()
{
  if (this->index!=NULL)
    delete this->index;
}

void ImportanceSampler::operator=(const ImportanceSampler &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->index!=NULL)
    delete this->index;

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* ImportanceSampler::smc_duplicate() const
{
  return( new ImportanceSampler(*this));
}

LikelihoodEstimator* ImportanceSampler::duplicate() const
{
  return( new ImportanceSampler(*this));
}

void ImportanceSampler::make_copy(const ImportanceSampler &another)
{
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

// void ImportanceSampler::smc_step()
// {
// }
//
// void ImportanceSampler::weight_update()
// {
// }


SMCOutput* ImportanceSampler::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  //std::cout<<simulation->back().back()->back().parameters<<std::endl;
  this->evaluate_smc(simulation);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* ImportanceSampler::initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void ImportanceSampler::simulate_smc(SMCOutput* current_state)
{
  this->simulate_proposal(current_state);
}

void ImportanceSampler::evaluate_smc(SMCOutput* current_state)
{
  //std::cout<<current_state->back().back()->back().parameters<<std::endl;
  //this->weight(current_state, conditioned_on_parameters);
  this->the_worker->weight(this->index,
                           current_state->back());
  //current_state->initialise_next_step();
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  //this->smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
}

void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  //this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_state->back());
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

MoveOutput* ImportanceSampler::move(RandomNumberGenerator &rng,
                                    Particle &particle)
{
  return new SinglePointMoveOutput(particle);
}

void ImportanceSampler::weight_for_adapting_sequence(Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles);
}

SMCOutput* ImportanceSampler::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* ImportanceSampler::initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void ImportanceSampler::simulate_smc(SMCOutput* current_state,
                                     const Parameters &conditioned_on_parameters)
{
  this->simulate_proposal(current_state, conditioned_on_parameters);
}

void ImportanceSampler::evaluate_smc(SMCOutput* current_state,
                                     const Parameters &conditioned_on_parameters)
{
  //this->weight(current_state, conditioned_on_parameters);
  this->the_worker->weight(this->index,
                           current_state->back(),
                           conditioned_on_parameters);
  //current_state->initialise_next_step();
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

void ImportanceSampler::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                const Parameters &conditioned_on_parameters)
{
  //this->smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back(),
                                    conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void ImportanceSampler::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                           const Parameters &conditioned_on_parameters)
{
  //this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_state->back(),
                                                      conditioned_on_parameters);
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

void ImportanceSampler::subsample_simulate_smc(SMCOutput* current_state,
                                     const Parameters &conditioned_on_parameters)
{
  this->simulate_proposal(current_state, conditioned_on_parameters);
}

void ImportanceSampler::subsample_evaluate_smc(SMCOutput* current_state,
                                               const Parameters &conditioned_on_parameters)
{
  //this->weight(current_state, conditioned_on_parameters);
  this->the_worker->subsample_weight(this->index,
                                     current_state->back(),
                                     conditioned_on_parameters);
  //current_state->initialise_next_step();
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

void ImportanceSampler::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                   const Parameters &conditioned_on_parameters)
{
  //this->smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->subsample_smcfixed_weight(this->index,
                                              current_state->back(),
                                              conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void ImportanceSampler::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                     const Parameters &conditioned_on_parameters)
{
  //this->smcadaptive_given_smcfixed_weight(current_state, conditioned_on_parameters);
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_state->back(),
                                                                conditioned_on_parameters);
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  
  //this->log_likelihood = current_state->log_likelihood_pre_last_step + current_state->latest_log_normalising_constant_ratio();
  
  current_state->log_likelihood = current_state->latest_log_normalising_constant_ratio();
}

MoveOutput* ImportanceSampler::move(RandomNumberGenerator &rng,
                                    Particle &particle,
                                    const Parameters &conditioned_on_parameters)
{
  return new SinglePointMoveOutput(particle);
}

void ImportanceSampler::weight_for_adapting_sequence(Particles &current_particles,
                                                     const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles,
                                                      conditioned_on_parameters);
}

MoveOutput* ImportanceSampler::subsample_move(RandomNumberGenerator &rng,
                                              Particle &particle,
                                              const Parameters &conditioned_on_parameters)
{
  return new SinglePointMoveOutput(particle);
}

void ImportanceSampler::subsample_weight_for_adapting_sequence(Particles &current_particles,
                                                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_particles,
                                                                conditioned_on_parameters);
}

/*
void ImportanceSampler::weight(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->weight(conditioned_on_parameters);
  current_state->initialise_next_step();
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
}

void ImportanceSampler::smcfixed_weight(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcfixed_weight(conditioned_on_parameters);
  current_state->initialise_next_step();
}


void ImportanceSampler::smcadaptive_given_smcfixed_weight(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
  current_state->weight_update(this->the_worker->get_unnormalised_log_incremental_weights());
}
 */

// Comment for later...
// IS: simulate
// SMC w MCMC: t=1 simulate, t>1 loop resample-move, then weight, until stopping point reached with resample-move being the final step
// PF/SMC: t=1 simulate, t>1 resample, sim prop, then weight, until stopping point reached with sim prop being final step
//void ImportanceSampler::smc_simulate(SMCOutput* current_state)
//{
   //this->output->add_proposed_particles(particles);

   //the_worker->simulate_and_weight();
  /*
   unsigned int number_of_points = algorithm["number_of_points"];

   List observed_data = model["observed_data"];

   // Do the initial importance sampling step.

   // Do the simulation.
   SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
   SimulateImportanceSamplingProposalPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

   std::vector<List> proposed_points;
   proposed_points.reserve(number_of_points);
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   proposed_points.push_back(simulate_proposal());
   }

   LikelihoodEstimator* likelihood_estimator = make_likelihood_estimator(model, algorithm);

   std::vector<List> proposed_auxiliary_variables;
   proposed_auxiliary_variables.reserve(number_of_points);

   for (std::vector<List>::const_iterator i=proposed_points.begin(); i!=proposed_points.end(); ++i)
   {
   proposed_auxiliary_variables.push_back(likelihood_estimator->simulate_auxiliary_variables(*i));
   }

   likelihood_estimator->is_setup_likelihood_estimator(proposed_points,
   proposed_auxiliary_variables);

   arma::colvec log_weights(number_of_points);
   bool prior_is_proposal = algorithm["prior_is_proposal"];
   if (prior_is_proposal==TRUE)
   {
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]);
   }
   }
   else
   {
   SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
   EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

   SEXP evaluate_log_proposal_SEXP = algorithm["evaluate_log_proposal"];
   EvaluateLogDistributionPtr evaluate_log_proposal = load_evaluate_log_distribution(evaluate_log_proposal_SEXP);

   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]) + evaluate_log_prior(proposed_points[i]) - evaluate_log_proposal(proposed_points[i]);
   }
   }

   if (likelihood_estimator != NULL)
   delete likelihood_estimator;

   return List::create(Named("proposed_points") = proposed_points,
   Named("proposed_auxiliary_variables") = wrap(proposed_auxiliary_variables),
   Named("log_weights") = log_weights,
   Named("log_normalising_constant") = log_sum_exp(log_weights));
   */

//}

//void ImportanceSampler::smc_weight(SMCOutput* current_state)
//{
//  this->the_worker->weight();
//  current_state->update_unnormalised_log_incremental_weights(this->get_unnormalised_log_incremental_weights());
  
//  current_state->update_unnormalised_log_weights();
//}
