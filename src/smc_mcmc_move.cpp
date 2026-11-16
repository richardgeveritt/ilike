#include <numeric>
#include <algorithm>
#include "smc_mcmc_move.h"
#include "chunk_parallel.h"
#include <mutex>
#include <memory>
#include "factors.h"
#include "factor_variables.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_index.h"

namespace ilike
{

// The factors in full_index that are not in used_index: those the SMC never
// evaluates because they cancel in its weights (typically the prior, when the
// prior is the initial proposal). Returns NULL if there are none. The caller
// owns the result.
static Index* make_cancelled_index(const Index* full_index,
                                   const Index* used_index)
{
  if ( (full_index==NULL) || (used_index==NULL) )
    return NULL;
  std::vector<size_t> cancelled;
  for (auto i=full_index->begin(); i!=full_index->end(); ++i)
  {
    if (std::find(used_index->begin(), used_index->end(), *i)==used_index->end())
      cancelled.push_back(*i);
  }
  if (cancelled.size()==0)
    return NULL;
  return new VectorIndex(cancelled);
}

SMCMCMCMove::SMCMCMCMove()
:SMC()
{
  this->mcmc = NULL;
  this->index = NULL;
  this->cancelled_index = NULL;
}

// Multiple MCMC chains.
SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         const std::vector<Parameters> &initial_points_in,
                         const arma::colvec &log_probabilities_of_initial_values_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     initial_points_in.size(),
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
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
  this->index = new VectorIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  this->particle_simulator = NULL;
  
  this->initial_particles = initial_points_in;
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = true;
}

// Multiple MCMC chains.
SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     number_of_particles_in,
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     true,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  this->proposed_particles_inputted = false;
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index = make_cancelled_index(full_index, without_cancelled_index);
  this->mcmc_at_last_step = true;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         const std::vector<Parameters> &initial_points_in,
                         const arma::colvec &log_probabilities_of_initial_values_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in,
     seed_in,
     data_in,
     algorithm_parameters,
     initial_points_in.size(),
     lag_in,
     lag_proposed_in,
     mcmc_in->get_proposals(),
     -1.0,
     false,
     true,
     true,
     transform_proposed_particles,
     results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  this->particle_simulator = NULL;
  
  this->initial_particles = initial_points_in;
  
  //arma::mat thing = this->initial_particles[0]["tau"];
  this->log_probabilities_of_initial_values = log_probabilities_of_initial_values_in;
  this->proposed_particles_inputted = true;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("SMCMCMCMove::SMCMCMCMove - need more than one initialising point.");
  
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              100,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index = make_cancelled_index(full_index, without_cancelled_index);
  this->mcmc_at_last_step = true;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), -1.0, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
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
  this->index = new VectorIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  SMCCriterion* smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

/*
 void SMCMCMCMove::set_multiple_mcmc(RandomNumberGenerator* rng_in,
 size_t* seed_in,
 const Data* data_in,
 MCMC* mcmc_in,
 SimulateDistributionPtr simulate_proposal_in,
 EvaluateLogDistributionPtr evaluate_log_prior_in,
 EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
 bool parallel_in,
 size_t grain_size_in)
 {
 this->model_and_algorithm.likelihood_estimators.resize(0);
 this->model_and_algorithm.likelihood_estimators.reserve(1);
 this->model_and_algorithm.likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
 seed_in,
 data_in,
 evaluate_log_prior_in,
 evaluate_log_likelihood_in,
 true));
 
 this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
 this->model_and_algorithm.likelihood_estimators);
 
 if (parallel_in==TRUE)
 {
 }
 else
 {
 this->the_worker = new SequentialSMCWorker(this,
 this->model_and_algorithm.particle_simulator,
 this->model_and_algorithm.likelihood_estimators);
 }
 
 std::vector<double> schedule_in;
 schedule_in.push_back(0.0);
 schedule_in.push_back(1.0);
 std::string variable_in = "power";
 this->model_and_algorithm.smc_criterion = new CESSSMCCriterion(-1.0);
 this->model_and_algorithm.smc_termination = NULL;
 this->sequencer = Sequencer(this->the_worker,
 schedule_in,
 variable_in,
 this->model_and_algorithm.smc_criterion,
 this->model_and_algorithm.smc_termination);
 
 this->mcmc = mcmc_in;
 }
 */

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::string variable_in = "power";
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  true);
  
  PowerFunctionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_times_likelihood,
                                                                  power,
                                                                  variable_in,
                                                                  false));
  indices.push_back(0);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   proposal,
                                                                                   true);
  
  PowerFunctionPtr second_power = annealing_one_minus_power;
  
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  variable_in,
                                                                  false));
  
  indices.push_back(1);
  this->index = new VectorIndex(indices);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         const std::vector<double> &temperatures_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateDistributionPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, true, true, true, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  std::string variable_in = "power";
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);
  //Parameters candidate_parameters = proposal->independent_simulate(*this->rng);
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  true);
  
  PowerFunctionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_times_likelihood,
                                                                  power,
                                                                  variable_in,
                                                                  false));
  indices.push_back(0);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   proposal,
                                                                                   true);
  
  PowerFunctionPtr second_power = annealing_one_minus_power;
  
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  variable_in,
                                                                  false));
  
  indices.push_back(1);
  this->index = new VectorIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  std::vector<double> schedule_in;
  schedule_in.reserve(temperatures_in.size()+1);
  schedule_in.push_back(0.0);
  for (size_t i=0; i<temperatures_in.size(); ++i)
  {
    schedule_in.push_back(temperatures_in[i]);
  }
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              number_of_bisections_in,
                              smc_criterion);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  VectorIndex mcmc_index(indices); // set_index copies it
  this->mcmc->set_index(&mcmc_index);
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         size_t number_of_bisections_in,
                         SMCTermination* termination_in,
                         const std::string &sequence_variable_in,
                         const std::vector<double> &schedule_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         bool proposal_is_evaluated_in,
                         bool smcfixed_flag_in,
                         bool sequencer_limit_is_fixed_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), resampling_desired_ess_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  //Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<size_t> indices;
  indices.reserve(likelihood_estimators_in.size());
  for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
    indices.push_back(i);
  this->index = new VectorIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  SMCCriterion* smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              sequence_variable_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  if (this->mcmc!=NULL)
  {
    VectorIndex mcmc_index(indices); // set_index_if_null copies it
    this->mcmc->set_index_if_null(&mcmc_index);
  }
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         const Parameters &algorithm_parameters,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         SMCCriterion* adaptive_resampling_in,
                         SMCCriterion* adaptive_target_in,
                         size_t number_of_bisections_in,
                         SMCTermination* termination_in,
                         const std::vector<std::string> &sequence_variables_in,
                         const std::vector<std::vector<double>> &schedules_in,
                         const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                         IndependentProposalKernel* proposal_in,
                         Index* without_cancelled_index,
                         Index* full_index,
                         bool proposal_is_evaluated_in,
                         bool smcfixed_flag_in,
                         bool sequencer_limit_is_fixed_in,
                         bool mcmc_at_last_step_in,
                         bool transform_proposed_particles,
                         bool parallel_in,
                         size_t grain_size_in,
                         const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, mcmc_in->get_proposals(), adaptive_resampling_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  this->chunk_size = grain_size_in;
  //Parameters candidate_parameters = proposal_in->independent_simulate(*this->rng);
  
  mcmc_in->set_proposal_parameters(&this->algorithm_parameters);
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
  //std::vector<size_t> indices;
  //indices.reserve(likelihood_estimators_in.size());
  //for (size_t i=0; i<likelihood_estimators_in.size(); ++i)
  //  indices.push_back(i);
  //this->index = new VectorIndex(indices);
  this->index = without_cancelled_index;
  
  this->factors = new VectorFactors(likelihood_estimators_in);
  
  /*
   for (auto i=likelihood_estimators.begin();
   i!=likelihood_estimators.end();
   ++i)
   {
   (*i)->setup(candidate_parameters);
   }
   */
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in);
  
  
  this->parallel = parallel_in;
  this->the_worker = new SequentialSMCWorker(this);
  
  SMCCriterion* smc_criterion = adaptive_target_in;
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedules_in,
                              sequence_variables_in,
                              number_of_bisections_in,
                              smc_criterion,
                              termination_in);
  //this->sequencer_parameters = &this->sequencer.schedule_parameters;
  
  this->mcmc = mcmc_in;
  this->cancelled_index = NULL;
  if (this->mcmc!=NULL)
    this->mcmc->set_index_if_null(full_index);
  this->cancelled_index = make_cancelled_index(full_index, without_cancelled_index);
  
  this->mcmc_at_last_step = mcmc_at_last_step_in;
}

//Copy constructor for the SMCMCMCMove class.
SMCMCMCMove::SMCMCMCMove(const SMCMCMCMove &another)
:SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCMCMCMove class.
SMCMCMCMove::~SMCMCMCMove()
{
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;
  
  if (this->cancelled_index!=NULL)
    delete this->cancelled_index;
}

void SMCMCMCMove::operator=(const SMCMCMCMove &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;
  
  if (this->cancelled_index!=NULL)
    delete this->cancelled_index;
  
  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCMCMCMove::smc_duplicate(void) const
{
  return( new SMCMCMCMove(*this));
}

LikelihoodEstimator* SMCMCMCMove::duplicate() const
{
  return( new SMCMCMCMove(*this));
}

void SMCMCMCMove::make_copy(const SMCMCMCMove &another)
{
  if (another.mcmc!=NULL)
    this->mcmc = another.mcmc->mcmc_duplicate();
  else
    this->mcmc = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
  
  if (another.cancelled_index!=NULL)
    this->cancelled_index = another.cancelled_index->duplicate();
  else
    this->cancelled_index = NULL;
  
  this->mcmc_at_last_step = another.mcmc_at_last_step;
}

SMCOutput* SMCMCMCMove::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

/*
 SMCOutput* SMCMCMCMove::specific_run(const std::string &directory_name)
 {
 SMCOutput* simulation = this->initialise_smc();
 this->simulate_smc(simulation);
 this->evaluate_smc(simulation);
 simulation->normalise_and_resample_weights();
 return simulation;
 }
 */

SMCOutput* SMCMCMCMove::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

int SMCMCMCMove::rich_chunk_mode() const
{
  if ((dynamic_cast<const SequentialSMCWorker*>(this->the_worker) == NULL) ||
      (this->mcmc == NULL) || (this->factors == NULL))
    return 0;
  if (this->factors->supports_parallel_chunks())
    return 1;
  if (this->factors->supports_parallel_copies())
    return 2;
  return 0;
}

void SMCMCMCMove::move_rich_in_chunks(Particles* next_particles,
                                      Particles* current_particles,
                                      bool copies)
{
  // SequentialSMCWorker::specific_move and the seed bookkeeping of
  // SMCWorker::move, grain_size particles at a time, each chunk with its own
  // stream, each particle's MoveOutput put in its own slot. A kernel that
  // adapts is copied for every chain (C98), which touches a static counter
  // of the kernels: such runs are chunked but not parallel.
  SMCWorker* worker = this->the_worker;
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  bool in_parallel = this->parallel && !this->mcmc->adapts();
  if (this->parallel && !in_parallel)
    this->parallel_refused = true;
  next_particles->particles.assign(N, NULL);

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  run_chunk_groups(0, K, in_parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     // the group's copy of the factors, drawing from its own
                     // generator and seed counter
                     std::unique_ptr<Factors> own_factors;
                     RandomNumberGenerator group_rng;
                     size_t group_seed = 0;
                     const Factors* factors_for_group = this->factors;
                     if (copies)
                     {
                       own_factors.reset(this->factors->duplicate());
                       own_factors->redirect_random_numbers(&group_rng, &group_seed);
                       factors_for_group = own_factors.get();
                     }
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                     {
                       size_t first = k*B;
                       size_t n = std::min(B, N - first);
                       RandomNumberGenerator local_rng(streams[k]);
                       for (size_t i = first; i < first + n; ++i)
                       {
                         const Particle &ancestor = (*current_particles)[current_particles->ancestor_variables[i]]->back();
                         MoveOutput* moved;
                         if (copies)
                         {
                           group_seed = size_t(local_rng());
                           Particle start = ancestor;
                           if ((start.factor_variables != NULL) && !start.factor_variables->repoint_to(factors_for_group))
                             Rcpp::stop("SMCMCMCMove::move_rich_in_chunks: the factors could not be repointed.");
                           moved = this->move(local_rng, start);
                           // what is kept outlives the group's copy
                           if (!moved->repoint_to(this->factors))
                             Rcpp::stop("SMCMCMCMove::move_rich_in_chunks: the factors could not be repointed.");
                         }
                         else
                         {
                           moved = this->move(local_rng, ancestor);
                         }
                         moved->back().previous_self = &ancestor;
                         next_particles->particles[i] = moved;
                       }
                     }
                   });

  // as SMCWorker::move
  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

bool SMCMCMCMove::can_move_flat() const
{
  // Conditions (dev/plan.md 2.6a, G3): the moves read the cached factor
  // values (so the cancelled factors are evaluated before each run, as in
  // move()); one kernel shared by all chains (no adaptation); a run length
  // fixed in advance (iterations(n): no termination reading the chain); and
  // every factor restorable from its summary (the SMC's factors) or
  // reinitialisable (the cancelled ones).
  if (this->force_rich_generations)
    return false;
  if (dynamic_cast<const SequentialSMCWorker*>(this->the_worker) == NULL)
    return false;
  if ((this->mcmc == NULL) || !this->mcmc->reads_cached_factors() || this->mcmc->adapts())
    return false;
  if (this->mcmc->known_number_of_iterations() == 0)
    return false;
  if ((this->factors == NULL) || !this->factors->supports_flat_moves(this->index, this->cancelled_index))
    return false;
  return true;
}

void SMCMCMCMove::restore_particle(Particle &particle,
                                   const Particles &generation,
                                   size_t row,
                                   const Parameters &restore_parameters) const
{
  // What a copy of the rich ancestor would be, for everything the moves
  // read: its variables (the SMC sequence variables stay shared with the
  // sequencer, as in every particle of this SMC), its factor outputs as its
  // last evaluation left them, and a new proposal store (gradients and
  // transforms are recomputed from the variables when needed).
  if (particle.parameters.vector_size() == 0)
    particle.parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  generation.flat_final_points.assign_into(row, particle.parameters);
  if (particle.factor_variables == NULL)
    particle.factor_variables = this->factors->simulate_factor_variables(particle.parameters);
  particle.factor_variables->set_particle(&particle);
  if (!particle.factor_variables->restore_from_summary(this->index,
                                                      this->cancelled_index,
                                                      generation.adaptive_summaries.colptr(row),
                                                      restore_parameters))
    Rcpp::stop("SMCMCMCMove::restore_particle: the factors of a particle could not be restored.");
  if (particle.ensemble_factor_variables != NULL)
  {
    delete particle.ensemble_factor_variables;
    particle.ensemble_factor_variables = NULL;
  }
  particle.current_proposal_store.clear();
  particle.simulate_proposal_variables(&this->proposals_to_transform_for,
                                       &this->proposals_to_find_gradient_for);
  particle.accepted_outputs.clear();
  particle.previous_self = NULL;
  particle.previous_target_evaluated = generation.adaptive_previous_targets[row];
  particle.target_evaluated = generation.adaptive_previous_targets[row];
  particle.subsample_previous_target_evaluated = 0.0;
  particle.subsample_target_evaluated = 0.0;
  particle.previous_ensemble_target_evaluated = 0.0;
  particle.subsample_previous_ensemble_target_evaluated = 0.0;
  particle.ensemble_target_evaluated = 0.0;
  particle.subsample_ensemble_target_evaluated = 0.0;
}

void SMCMCMCMove::move_flat(Particles* next_particles,
                            Particles* current_particles)
{
  // The steps of SequentialSMCWorker::specific_move with this->move() and
  // MCMC::run for each particle, then the seed bookkeeping of SMCWorker::move.
  // The chunks run in order, or in parallel (G4, G5), each with its own
  // stream and writing only its own particles' rows: a chain of
  // `iterations` moves stores iterations/thin states, so particle i's states
  // are rows [i*m, (i+1)*m) of the chain store.
  SMCWorker* worker = this->the_worker;
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  size_t iterations = this->mcmc->known_number_of_iterations();
  size_t thin = std::max<size_t>(this->mcmc->get_thin(), 1);
  size_t m = iterations/thin;
  bool from_flat = current_particles->is_flat();
  // the factor values a rich generation's particles held after the target
  // was chosen: those of the last trial target (SMC::finish_weighting_from_summaries)
  const Parameters &restore_parameters =
    (from_flat && current_particles->adaptive_summaries_used) ? current_particles->last_trial_schedule
                                                               : this->sequencer.schedule_parameters;

  next_particles->flat = true;
  next_particles->flat_size = N;
  next_particles->flat_output_lengths.assign(N, double(m));
  next_particles->adaptive_previous_targets.set_size(N);

  bool keep_records = this->write_factor_output;
  size_t number_of_records = 0;
  int summary_size = -1;
  // the first state to be stored, and the first particle to finish, fix the
  // layouts (every particle gives the same; set_point checks it)
  std::once_flag chain_setup;
  std::once_flag final_setup;

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  auto run_chunk = [&](size_t k, Particle &current, Particle &scratch, std::vector<double> &records)
  {
    size_t first = k*B;
    size_t n = std::min(B, N - first);
    RandomNumberGenerator local_rng(streams[k]);

    for (size_t i = first; i < first + n; ++i)
    {
      size_t ancestor = current_particles->ancestor_variables[i];
      if (from_flat)
        this->restore_particle(current, *current_particles, ancestor, restore_parameters);
      else
        current = (*current_particles)[ancestor]->back();

      // as move(): the cancelled factors, then the run
      if (this->cancelled_index != NULL)
        current.evaluate_likelihoods(this->cancelled_index);

      size_t stored = 0;
      for (size_t iteration = 1; iteration <= iterations; ++iteration)
      {
        this->mcmc->move_in_place(local_rng, current, scratch);
        current.tell_factors_to_forget_they_were_already_written_to_file();
        if ((thin > 1) && ((iteration % thin) != 0))
          continue;
        // as StandardMCMCOutput::record
        std::call_once(chain_setup, [&]()
        {
          if (keep_records && (current.factor_variables != NULL))
          {
            if (!current.factor_variables->has_compact_values())
              Rcpp::stop("SMCMCMCMove::move_flat: a factor has no compact record.");
            number_of_records = current.factor_variables->number_of_compact_values();
            next_particles->flat_points.set_number_of_extra_columns(number_of_records);
            next_particles->flat_record_writer = current.factor_variables->duplicate();
            next_particles->flat_record_writer->set_particle(NULL);
          }
          next_particles->flat_points.prepare(current.parameters, N*m);
        });
        size_t row = i*m + stored;
        if (number_of_records > 0)
        {
          records.resize(number_of_records);
          current.factor_variables->compact_values(records.data());
          next_particles->flat_points.set_point(row, current.parameters, records.data());
        }
        else
        {
          next_particles->flat_points.set_point(row, current.parameters);
        }
        ++stored;
      }

      // as move(): the target over the SMC's index, which becomes the
      // particle's previous target (set_previous_target_evaluated_to_target_evaluated)
      double target = current.get_cached_likelihoods(this->index);
      const FactorVariables* factor_variables = current.factor_variables;
      int size = (factor_variables == NULL) ? -1 : factor_variables->adaptive_summary_size(this->index);
      std::call_once(final_setup, [&]()
      {
        next_particles->flat_final_points.prepare(current.parameters, N);
        if (size < 0)
          Rcpp::stop("SMCMCMCMove::move_flat: a factor has no adaptive summary.");
        summary_size = size;
        next_particles->adaptive_summaries.set_size(size_t(summary_size), N);
        next_particles->flat_summary_evaluator = factor_variables->duplicate();
        next_particles->flat_summary_evaluator->set_particle(NULL);
        next_particles->adaptive_summary_evaluator = next_particles->flat_summary_evaluator;
        next_particles->adaptive_summary_index = this->index;
        next_particles->adaptive_summary_state = 1;
      });
      next_particles->flat_final_points.set_point(i, current.parameters);
      if (size != summary_size)
      {
        Rcpp::stop("SMCMCMCMove::move_flat: the particles' adaptive summaries differ in size.");
      }
      if (summary_size > 0)
        factor_variables->adaptive_summary(this->index, next_particles->adaptive_summaries.colptr(i));
      next_particles->adaptive_previous_targets[i] = target;
    }
  };


  // the moves read the kernels and factors only (C98; checked by
  // supports_flat_moves, which admits exact and annealed-exact factors only)
  bool in_parallel = this->parallel && (this->factors != NULL) && this->factors->supports_parallel_chunks();
  if (this->parallel && !in_parallel)
    this->parallel_refused = true;
  run_chunk_groups(0, K, in_parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     Particle current;
                     Particle scratch;
                     std::vector<double> records;
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                       run_chunk(k, current, scratch, records);
                   });

  // as SMCWorker::move
  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    if (this->can_move_flat())
    {
      this->move_flat(next_particles,
                      current_particles);
    }
    else if (this->rich_chunk_mode() > 0)
    {
      this->move_rich_in_chunks(next_particles,
                                current_particles,
                                this->rich_chunk_mode() == 2);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->the_worker->move(next_particles,
                             current_particles);
    }
    
    //this->evaluate_smcfixed_part_smc(current_state);
    current_state->increment_smc_iteration();
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state)
{
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    Rcpp::stop("SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,this->index);
    this->finish_weighting_from_summaries(this->index, current_state->back());
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    // One log-likelihood and one time per stored particle set: recorded here,
    // before the termination check, so that the particle set written out by
    // the final MCMC move below has both.
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time_and_reset_start();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state);
        current_state->set_llhd(current_state->log_likelihood);
        current_state->set_time_and_reset_start();
        //current_state->decrement_smc_iteration();
      }
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state);
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  // The MCMC moves read their acceptance-ratio denominators from the
  // particle's per-factor cache, so every factor in the target must have been
  // evaluated at the current parameters. The SMC has evaluated the factors in
  // its own index; evaluate the ones it skipped (cancelled_index) here. Only
  // those are evaluated, so random (pseudo-marginal) factors the SMC has
  // already estimated are not re-estimated.
  // A chain of unadjusted moves only never reads cached values, and the user
  // may not have supplied the functions to evaluate the factors, so nothing is
  // evaluated here (the behaviour before the cache was used).
  if (!this->mcmc->reads_cached_factors())
    return this->mcmc->run(rng, particle);
  
  MoveOutput* output;
  if (this->cancelled_index!=NULL)
  {
    Particle start = particle;
    start.evaluate_likelihoods(this->cancelled_index);
    output = this->mcmc->run(rng, start);
  }
  else
  {
    output = this->mcmc->run(rng, particle);
  }
  
  // The SMC uses target_evaluated (over its own index) in its next weight
  // update, so set it from the cache rather than leaving whatever the last
  // move happened to store there.
  output->back().target_evaluated = output->back().get_cached_likelihoods(this->index);
  return output;
}

//void SMCMCMCMove::weight_for_adapting_sequence(Particles &current_particles)
//{
//  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
//                                                      current_particles);
//}

SMCOutput* SMCMCMCMove::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}

/*
 SMCOutput* SMCMCMCMove::specific_run(const std::string &directory_name,
 const Parameters &conditioned_on_parameters)
 {
 SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
 this->simulate_smc(simulation, conditioned_on_parameters);
 this->evaluate_smc(simulation, conditioned_on_parameters);
 simulation->normalise_and_resample_weights();
 return simulation;
 }
 */

SMCOutput* SMCMCMCMove::specific_initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->move(next_particles,
                           current_particles);
    
    //this->evaluate_smcfixed_part_smc(current_state,
    //                                 conditioned_on_parameters);
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  /*
   if (this->sequencer_parameters!=NULL)
   {
   Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
   this->the_worker->smcfixed_weight(this->index,
   current_state->back(),
   all_parameters);
   }
   else
   {
   this->the_worker->smcfixed_weight(this->index,
   current_state->back(),
   conditioned_on_parameters);
   }
   */
  
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                               const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index);
    this->finish_weighting_from_summaries(this->index, current_state->back());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_time_and_reset_start();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        this->simulate_smc(current_state,conditioned_on_parameters);
        current_state->set_llhd(current_state->log_likelihood);
        current_state->set_time_and_reset_start();
        //current_state->decrement_smc_iteration();
      }
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::subsample_simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    //this->subsample_evaluate_smcfixed_part_smc(current_state);
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::subsample_simulate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    //this->subsample_evaluate_smcfixed_part_smc(current_state,
    //                                           conditioned_on_parameters);
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::subsample_evaluate_smc(SMCOutput* current_state)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMCMCMove::subsample_evaluate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  /*
   if (this->sequencer_parameters!=NULL)
   {
   Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
   this->the_worker->subsample_smcfixed_weight(this->index,
   current_state->back(),
   all_parameters);
   }
   else
   {
   this->the_worker->subsample_smcfixed_weight(this->index,
   current_state->back(),
   conditioned_on_parameters);
   }
   */
  this->the_worker->subsample_smcfixed_weight(this->index,
                                              current_state->back());
  //current_state->initialise_next_step();
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
  /*
   if (this->sequencer_parameters!=NULL)
   {
   Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
   this->the_worker->subsample_smcfixed_weight(this->index,
   current_state->back(),
   all_parameters);
   }
   else
   {
   this->the_worker->subsample_smcfixed_weight(this->index,
   current_state->back(),
   conditioned_on_parameters);
   }
   */
  this->the_worker->subsample_smcfixed_weight(this->index,
                                              current_state->back());
  //current_state->initialise_next_step();
}

void SMCMCMCMove::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  
  // need to do something with sequencer for each particle
  
  // either set sequencer with paraneters - for now, just setting an end point is fine, but could make it more sophisticated
  
  // or reset sequencer
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index);
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        // The final move targets the same (subsampled) distribution as the
        // other moves of this run.
        this->subsample_simulate_smc(current_state);
        current_state->set_llhd(current_state->subsample_log_likelihood);
        current_state->set_time_and_reset_start();
        //current_state->decrement_smc_iteration();
      }
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                         const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // need to do something with sequencer for each particle
  
  // either set sequencer with paraneters - for now, just setting an end point is fine, but could make it more sophisticated
  
  // or reset sequencer
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index);
    
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      if (this->mcmc_at_last_step)
      {
        // The final move targets the same (subsampled) distribution as the
        // other moves of this run.
        this->subsample_simulate_smc(current_state, conditioned_on_parameters);
        current_state->set_llhd(current_state->subsample_log_likelihood);
        current_state->set_time_and_reset_start();
        //current_state->decrement_smc_iteration();
      }
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state,
                                 conditioned_on_parameters);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

/*
 MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 return this->mcmc->run(rng,
 particle,
 conditioned_on_parameters);
 }
 */

void SMCMCMCMove::weight_for_adapting_sequence(const Index* index,
                                               Particles &current_particles)
{
  // from the generation's adaptive summaries when available (dev/plan.md
  // 2.6a, G2; the bisection calls this once per trial target)
  this->weight_from_summaries_if_possible(this->index,
                                          current_particles);
}

/*
 void SMCMCMCMove::weight_for_adapting_sequence(const Index* index,
 Particles &current_particles,
 const Parameters &conditioned_on_parameters)
 {
 this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
 current_particles,
 conditioned_on_parameters);
 }
 */

MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  return this->mcmc->subsample_run(rng,
                                   particle);
}

/*
 MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 return this->mcmc->subsample_run(rng,
 particle,
 conditioned_on_parameters);
 }
 */

void SMCMCMCMove::subsample_weight_for_adapting_sequence(const Index* index,
                                                         Particles &current_particles)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_particles);
}

/*
 void SMCMCMCMove::subsample_weight_for_adapting_sequence(const Index* index,
 Particles &current_particles,
 const Parameters &conditioned_on_parameters)
 {
 this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
 current_particles,
 conditioned_on_parameters);
 }
 */

// void SMCMCMCMove::smc_step(void)
// {
//
// }
//
// void SMCMCMCMove::weight_update(void)
// {
//
// }

//void SMCMCMCMove::smc_update(SMCOutput* current_state)
//{
/*
 unsigned int number_of_points = algorithm["number_of_points"];
 
 List observed_data = model["observed_data"];
 
 // Do the initial importance sampling step.
 
 // Do the simulation.
 SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
 SimulateDistributionPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);
 
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
}
